#ifndef OPM_H__
#define OPM_H__


#include <emp-tool/emp-tool.h>
#include "emp-ot/emp-ot.h"
#include <iostream>
#include <algorithm>
#include <random>

using namespace emp;
using std::future;

template<typename IO>
class OPM {
public: 
    int party, threads_n;
    int height, dimension;

    IO *io;

    // for OPM receiver
    std::vector<uint32_t> permutation;

    // for OPM sender
    block *ggm_master_key = nullptr; 

    block *opm_sen_memery = nullptr, *opm_rev_memory = nullptr; //length = dimension * (2*dimension). //FIXME: will remove opm_rev_memory later
    block *column_sum; //length = dimension
    block *opm_a = nullptr, *opm_b = nullptr, *opm_delta = nullptr;
    
    block *ot_msg = nullptr, *ot_rev_msg = nullptr; //FIXME: will remove ot_rev_msg later

    PRG prg;


    OPM(int party, int threads_n, int height) {
        this->party = party;
        this->threads_n = threads_n;
        this->height = height;
        this->dimension = 1<<(height-1); // the last layer is used to check

        this->column_sum = new block[dimension]; 

        if (party == ALICE) receiver_init(); 
        else sender_init();
    }

    void receiver_init() {
        this->opm_rev_memory = new block[this->dimension * 2 * this->dimension]; 

        this->opm_delta = new block[this->dimension];
        for (int i = 0;  i < this->dimension; i++) {
            *(this->opm_delta + i) = zero_block;
        }

        this->ot_rev_msg = new block[this->height-1]; 
    }

    void sender_init() {
        this->ggm_master_key = new block[this->dimension]; 
        prg.random_block(this->ggm_master_key, this->dimension);

        this->opm_sen_memery = new block[this->dimension * 2 * this->dimension]; 

        this->opm_a = new block[this->dimension]; 
        this->opm_b = new block[this->dimension]; 
        for (int i = 0;  i < this->dimension; i++) {
            *(this->opm_a + i) = zero_block;
            *(this->opm_b + i) = zero_block;
        }
        
        this->ot_msg = new block[(this->height - 1)* 2];
    }

    ~OPM() {
        if(ggm_master_key != nullptr) delete[] ggm_master_key;

        delete[] column_sum;

        if (opm_sen_memery != nullptr) delete[] opm_sen_memery;
        if (opm_rev_memory != nullptr) delete[] opm_rev_memory;

        if(opm_a!=nullptr) delete[] opm_a;
        if(opm_b!=nullptr) delete[] opm_b;
        if(opm_delta!=nullptr) delete[] opm_delta;
        if(ot_msg != nullptr) delete[] ot_msg;
        if(ot_rev_msg != nullptr) delete[] ot_rev_msg; 
    }


    //TODO: use cryptigraphic secure random source
    void random_permutation() {
        for (int i = 0; i < this->dimension; i++) 
            this->permutation.push_back(i);
        
        std::random_device rand_div("/dev/urandom");
        auto rng = std::default_random_engine {rand_div()};
        std::shuffle(permutation.begin(), permutation.end(), rng);
    }

    void OPM_gen() { //`ot_msg_0` and `ot_msg_1` is of length `n * (height - 1)`

        block *ot_msg_0 = this->ot_msg, *ot_msg_1 = this->ot_msg + this->height - 1;  

        for (int i = 0; i < this->dimension; i++) {
            ggm_tree_gen(this->ggm_master_key + i, this->height, ot_msg_0 + i*(this->height-1), ot_msg_1 + i*this->height, this->opm_sen_memery + i*(2*this->dimension));
        }

        for (int i = 0; i < this->dimension; i++) {
            for(int j = 0; j < this->dimension; j++) {
                *(this->column_sum + j) = *(this->column_sum + j) ^ (*(this->opm_sen_memery + 2*i*this->dimension + 2*j));
                *(this->opm_b + i) = *(this->opm_b + i) ^ (*(this->opm_sen_memery + 2*i*this->dimension + 2*j + 1));
                *(this->opm_a + j) = *(this->opm_a + j) ^ (*(this->opm_sen_memery + 2*i*this->dimension + 2*j + 1));
            }
        }
    }

    void ggm_tree_gen(block seed, int depth, block *ot_msg_0, block *ot_msg_1, block *ggm_tree) {
		//this->ggm_tree = ggm_tree_mem; // this->ggm_tree points to `ggm_tree_mem` from the outter caller
		TwoKeyPRP *prp = new TwoKeyPRP(zero_block, makeBlock(0, 1));
		prp->node_expand_1to2(ggm_tree, seed); // expand `seed` to get two children and store at ggm_tree[0] and ggm_tree[1]
		ot_msg_0[0] = ggm_tree[0];
		ot_msg_1[0] = ggm_tree[1];
		prp->node_expand_2to4(&ggm_tree[0], &ggm_tree[0]); /* expand ggm_tree[0] and ggm_tree[1] respectively, then store four extended 
															  children sequentically, starting from ggm_tree[0]. 
															  *ALL CHILDREN NODES ARE STORED IN THE FRONT*.
														    */
		ot_msg_0[1] = ggm_tree[0] ^ ggm_tree[2];
		ot_msg_1[1] = ggm_tree[1] ^ ggm_tree[3];
		for(int h = 2; h < depth-1; ++h) {
			ot_msg_0[h] = ot_msg_1[h] = zero_block;
			int sz = 1<<h;
			for(int i = sz-4; i >=0; i-=4) {
				prp->node_expand_4to8(&ggm_tree[i*2], &ggm_tree[i]);
				ot_msg_0[h] ^= ggm_tree[i*2];
				ot_msg_0[h] ^= ggm_tree[i*2+2];
				ot_msg_0[h] ^= ggm_tree[i*2+4];
				ot_msg_0[h] ^= ggm_tree[i*2+6];
				ot_msg_1[h] ^= ggm_tree[i*2+1];
				ot_msg_1[h] ^= ggm_tree[i*2+3];
				ot_msg_1[h] ^= ggm_tree[i*2+5];
				ot_msg_1[h] ^= ggm_tree[i*2+7];
			}
		}
		delete prp;
	}

    void ggm_tree_reconstruction(bool *b, block *m, int depth, block *ggm_tree) {
		int to_fill_idx = 0;
		TwoKeyPRP prp(zero_block, makeBlock(0, 1));
		for(int i = 1; i < depth; ++i) {
			to_fill_idx = to_fill_idx * 2;
			ggm_tree[to_fill_idx] = ggm_tree[to_fill_idx+1] = zero_block; // two undefined children nodes, one will be filled
			if(b[i-1] == false) {// `b[i-1] = 0` <=> `choice_pos[i] = 1`
				layer_recover(i, 0, to_fill_idx, m[i-1], &prp, ggm_tree);
				to_fill_idx += 1;
			} else layer_recover(i, 1, to_fill_idx+1, m[i-1], &prp, ggm_tree);
		}
	}

	void layer_recover(int depth, int lr, int to_fill_idx, block sum, TwoKeyPRP *prp, block *ggm_tree) {
		int layer_start = 0;
		int item_n = 1<<depth;
		block nodes_sum = zero_block;
		int lr_start = lr==0?layer_start:(layer_start+1);
		
		for(int i = lr_start; i < item_n; i+=2)
			nodes_sum = nodes_sum ^ ggm_tree[i];
		ggm_tree[to_fill_idx] = nodes_sum ^ sum;
		if(depth == this->height-1) return;
		if(item_n == 2)
			prp->node_expand_2to4(&ggm_tree[0], &ggm_tree[0]);
		else {
			for(int i = item_n-4; i >= 0; i-=4)
				prp->node_expand_4to8(&ggm_tree[i*2], &ggm_tree[i]);
		}
	}

    void get_choice_bits(int position, bool *b, int length) { // length >= 2
		for(int i = 0; i < length-1; i++) {
			b[length-2-i] = 1 - position % 2;
            position >>= 1;
		}
	}

    void OPM_recover(block *ot_rev_msg) { //recover OPM from punctured keys
        bool *b = new bool[2*this->dimension]; 
        for (int i = 0; i < this->dimension; i++) { //for the i-th ggm tree
            get_bits(2*this->permutation[i], b, this->height);
            ggm_tree_reconstruction(b, this->ot_rev_msg + i*(this->height-1) , this->height, this->opm_rev_memory + i*2*this->dimension);
        }
        delete[] b;
    }
    
    void test_ggm_tree(int depth, int index) {
        int n = 1<<(depth-1);
        block seed;
        block ggm_tree_sen[n];
        block ggm_tree_rev[n];
        block ot_msg_sen[2*(depth-1)];
        block ot_msg_rev[depth-1];
        bool b[depth-1];

        prg.random_block(&seed, 1);
        ggm_tree_gen(seed, depth, ot_msg_sen, ot_msg_sen+depth-1, ggm_tree_sen);
        
        std::cout<<"\033[32m"<< "before get_bits." <<std::endl;
        get_choice_bits(index, b, depth);
        for(int i = 0; i < depth-1; i++) {//simulating OTs
            std::cout<< b[i] << " ";
            ot_msg_rev[i] = (b[i] == false ? ot_msg_sen[i] : ot_msg_sen[depth-1+i]);
        }
        std::cout<<std::endl;

        ggm_tree_reconstruction(b, ot_msg_rev, depth, ggm_tree_rev);
        std::cout<< "n = " << n << std::endl;
        for(int i = 0; i < n; i++) {
            if(i != index && memcmp(&ggm_tree_sen[i], &ggm_tree_rev[i], sizeof(block)) != 0)
            {
                std::cout<< "ERROR - "<< i << " - " << ggm_tree_sen[i] <<" : "<< ggm_tree_rev[i] <<std::endl;
            }
        }

        std::cout<<"\033[32m"<< "before delete." <<std::endl;

        // delete[] ggm_tree_sen;
        // delete[] ggm_tree_rev;
        // delete[] ot_msg_sen;
        // delete[] ot_msg_rev;
        // delete[] b;
    }
};

#endif 