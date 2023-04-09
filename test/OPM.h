#ifndef OPM_H__
#define OPM_H__


#include <emp-tool/emp-tool.h>
#include "emp-ot/emp-ot.h"
#include "test/shuffle_common.h"
#include <iostream>
#include <algorithm>
#include <random>

using namespace emp;
using std::future;

//#define MALICIOUS


template<typename IO>
class OPM {
public: 
    int party, threads_n;
    size_t height, dimension;

    IO *io;
    OT<IO> *ot;

    // for OPM receiver
    std::vector<size_t> permutation, permutation_inverse;

    // for OPM sender
    block *ggm_master_key = nullptr; 

    block *opm_sen_memory = nullptr, *opm_rev_memory = nullptr; // length = dimension * (2*dimension).
    block *column_sum, *opm_tag;
    block *opm_a = nullptr, *opm_b = nullptr, *opm_delta = nullptr;
    
    block *ot_sen_msg = nullptr, *ot_rev_msg = nullptr;
    bool  *ot_choise_bits = nullptr;

    PRG prg;


<<<<<<< HEAD
    OPM(int party, IO *io, OT<IO> *ot, int threads_n, int height, bool pre_setup = false) { // height > 3 for OPM
=======
    OPM(int party, int threads_n, int height) { // height > 3 for OPM
>>>>>>> f9de09bd64820fe22cff2ec1e66b1769ab811ef0
        this->party = party;
        this->io = io;
        this->ot = ot;
        this->threads_n = threads_n;
        this->height = height;
        this->dimension = 1<<(height-2); // the last layer is used for checking, thus only layers `dimension= #nodes of layer (height-1)`, which is `1 <<(height-2)`.  
        this->column_sum = new block[dimension];
<<<<<<< HEAD
        memset(this->column_sum, 0, this->dimension*sizeof(block)); 
        this->opm_tag = new block;
        // if (party == BOB) receiver_init(); 
        // else sender_init();

        receiver_init(pre_setup);
        sender_init(pre_setup);
    }

    void receiver_init(bool pre_setup = false) {
        // std::cout<< "new opm_rev_memory" <<std::endl;
        this->opm_rev_memory = new block[2 * this->dimension * this->dimension]; 
        memset(this->opm_rev_memory, 0, 2 * this->dimension * this->dimension * sizeof(block));

        this->opm_delta = new block[this->dimension];
        memset(this->opm_delta, 0, this->dimension*sizeof(block)); 

        ot_rev_msg = new block[(height-2)*dimension];
        ot_choise_bits = new bool[(height-2)*dimension];
        
        if(pre_setup) { // setup receiver from local stored file

            return;
        }

        this->random_permutation();
    }

    void sender_init(bool pre_setup = false) {
=======
        this->opm_tag = new block;
        // if (party == ALICE) receiver_init(); 
        // else sender_init();
        // FIXME: will do initilization depending on party's role
        receiver_init();
        sender_init();
    }

    void receiver_init() {
        // std::cout<< "new opm_rev_memory" <<std::endl;
        this->opm_rev_memory = new block[2 * this->dimension * this->dimension]; 

        this->opm_delta = new block[this->dimension];
        memset(this->opm_delta, 0, this->dimension*sizeof(block)); 

        ot_rev_msg = new block[(height-2)*dimension];
        ot_choise_bits = new bool[(height-2)*dimension];
        
        this->random_permutation(); 
    }

    void sender_init() {
>>>>>>> f9de09bd64820fe22cff2ec1e66b1769ab811ef0

        this->ggm_master_key = new block[this->dimension];
        prg.random_block(this->ggm_master_key, this->dimension);

        this->opm_sen_memory = new block[2 * this->dimension * this->dimension]; 
        this->opm_a = new block[this->dimension]; 
        this->opm_b = new block[this->dimension];
<<<<<<< HEAD
        memset(this->opm_sen_memory, 0, 2 * this->dimension * this->dimension * sizeof(block));
        memset(this->opm_a, 0, this->dimension*sizeof(block));
        memset(this->opm_b, 0, this->dimension*sizeof(block));
        this->ot_sen_msg = new block[2*(this->height-2)*this->dimension];// layers [2, height-1] need OT messages
        if(pre_setup) { // setup sender from local stored file
            
            return;
        }
=======
        memset(this->opm_a, 0, this->dimension*sizeof(block));
        memset(this->opm_b, 0, this->dimension*sizeof(block));
        
        this->ot_sen_msg = new block[2*(this->height-2)*this->dimension];// layers [2, height-1] need OT messages
>>>>>>> f9de09bd64820fe22cff2ec1e66b1769ab811ef0
    }

    ~OPM() {
        if(ggm_master_key != nullptr) delete[] ggm_master_key;

        delete[] column_sum;
        delete opm_tag; 

        if (opm_sen_memory != nullptr) delete[] opm_sen_memory; 
        if (opm_rev_memory != nullptr) delete[] opm_rev_memory;

        if(opm_a != nullptr) delete[] opm_a;
        if(opm_b != nullptr) delete[] opm_b;
        if(opm_delta != nullptr) delete[] opm_delta;
        if(ot_sen_msg != nullptr) delete[] ot_sen_msg;
        if(ot_rev_msg != nullptr) delete[] ot_rev_msg;
        if(ot_choise_bits != nullptr) delete[] ot_choise_bits; 

<<<<<<< HEAD
        // std::cout<<"~OPM() end."<<std::endl;
=======
        std::cout<<"~OPM() end."<<std::endl;
>>>>>>> f9de09bd64820fe22cff2ec1e66b1769ab811ef0
    }


    //TODO: use cryptigraphic secure random source
    void random_permutation() {
        for (size_t i = 0; i < this->dimension; i++) { 
            this->permutation.push_back(i);
        }
        
        std::random_device rand_div("/dev/urandom");
        auto rng = std::default_random_engine {rand_div()};
        std::shuffle(permutation.begin(), permutation.end(), rng);

        permutation_inverse.resize(dimension);
        for(size_t i = 0; i< dimension; i++) {
            permutation_inverse[permutation[i]] = i;   
        }
    }

    void opm_gen() { //`ot_msg_0` and `ot_msg_1` is with length `n*(height-1)`
<<<<<<< HEAD

        block *ot_sen_msg_0 = this->ot_sen_msg, *ot_sen_msg_1 = this->ot_sen_msg + (this->height-2)*(this->dimension);

        // std::cout<< "   Before opv_gen()" <<std::endl;
        for (size_t i = 0; i < this->dimension; i++) { // TODO: parallel this 
            opv_gen(this->ggm_master_key[i], this->height, ot_sen_msg_0+i*(this->height-2), ot_sen_msg_1+i*(this->height-2), this->opm_sen_memory+2*i*this->dimension );
        }
        // std::cout<< "   After opv_gen()" <<std::endl;

        // TODO: parallel this
        memset(this->column_sum, 0, this->dimension*sizeof(block)); 
        for (size_t i = 0; i < this->dimension; i++) {
            for(size_t j = 0; j < this->dimension; j++) {
                this->column_sum[j] = this->column_sum[j] ^ this->opm_sen_memory[2*i*this->dimension + 2*j];
                this->opm_b[i] = this->opm_b[i] ^ this->opm_sen_memory[2*i*this->dimension + 2*j + 1]; 
                this->opm_a[j] = this->opm_a[j] ^ this->opm_sen_memory[2*i*this->dimension + 2*j + 1];
            }
        }
    }

    void opv_gen(block seed, int height, block *ot_msg_0, block *ot_msg_1, block *opv_memory) {
        //std::cout<< "in opv_gen() 1" <<std::endl;
		TwoKeyPRP *prp = new TwoKeyPRP(zero_block, makeBlock(0, 1));
		prp->node_expand_1to2(opv_memory, seed); //expand `seed` to get two children and store at ggm_tree[0] and ggm_tree[1]
        //std::cout<< "in opv_gen() 2" <<std::endl;
		ot_msg_0[0] = opv_memory[0];
		ot_msg_1[0] = opv_memory[1];
		prp->node_expand_2to4(&opv_memory[0], &opv_memory[0]); /* expand ggm_tree[0] and ggm_tree[1] respectively, then store four extended 
															  children sequentically, starting from ggm_tree[0]. 
															  *ALL CHILDREN NODES ARE STORED IN THE FRONT*.
														    */
        if((height-2) != 1){
            ot_msg_0[1] = opv_memory[0] ^ opv_memory[2];
		    ot_msg_1[1] = opv_memory[1] ^ opv_memory[3];
        }
		for(size_t h = 2; h < height-1; ++h) {
            if(h != (height-2)) {
			    ot_msg_0[h] = ot_msg_1[h] = zero_block;
            }
			size_t sz = 1<<h;
            //std::cout<< "\n\nh = " << h << ", sz = " << sz <<std::endl; 
			for(size_t i = sz-4; i >=0; i-=4) {
				prp->node_expand_4to8(&opv_memory[i*2], &opv_memory[i]);
                if(h != (height-2)) { //no OT messages needed for the last layer, column-wise sum will be sent for the last layer.
                    ot_msg_0[h] ^= opv_memory[i*2];
                    ot_msg_0[h] ^= opv_memory[i*2+2];
                    ot_msg_0[h] ^= opv_memory[i*2+4];
                    ot_msg_0[h] ^= opv_memory[i*2+6];
                    ot_msg_1[h] ^= opv_memory[i*2+1];
                    ot_msg_1[h] ^= opv_memory[i*2+3];
                    ot_msg_1[h] ^= opv_memory[i*2+5];
                    ot_msg_1[h] ^= opv_memory[i*2+7];
                }
                //std::cout<<"opv_memory[i*2+7] = " << i*2+7 << std::endl;
			}
		}
		delete prp;
	}

    void opm_reconstruction() { //`ot_msg_0` and `ot_msg_1` is with length `n*(height-1)`

        bool b[this->height-2];
        for(size_t i = 0; i < this->dimension; i++) { // TODO: parallel this 
            get_choice_bits(this->permutation[i], b, this->height-2);
            opv_reconstruction(b, this->ot_rev_msg+i*(this->height-2), this->height, opm_rev_memory+i*2*dimension);
            // log::green("\n\nnew recon:\n");
            // for(size_t j = 0; j < dimension; j++) {
            //     std::cout<<"opm_sen_memory[2*i*dimension+2*j]   = "<< opm_sen_memory[2*i*dimension+2*j] <<  ", opm_rev_memory[2*i*dimension+2*j]   = "<< opm_rev_memory[2*i*dimension+2*j] <<std::endl;
            //     std::cout<<"opm_sen_memory[2*i*dimension+2*j+1] = "<< opm_sen_memory[2*i*dimension+2*j+1] <<", opm_rev_memory[2*i*dimension+2*j+1] = "<< opm_rev_memory[2*i*dimension+2*j+1] <<std::endl;
            // }
            // std::cout<<std::endl;
        }
        
        //update opm_delta and prepare for opm sacrifice TODO: parallel this
        block *tmp_sum = new block[dimension]; 
        memset(tmp_sum, 0, sizeof(block)*dimension);
        for (size_t i = 0; i < dimension; i++) {
            for(size_t j = 0; j < dimension; j++) {
                tmp_sum[j] = tmp_sum[j] ^ opm_rev_memory[2*i*dimension + 2*j];
                if(j != permutation[i]) { 
                    opm_delta[i] ^= opm_rev_memory[2*i*dimension + 2*j+1];
                    opm_delta[permutation_inverse[j]] ^= opm_rev_memory[2*i*dimension + 2*j+1];
                }
=======

        block *ot_sen_msg_0 = this->ot_sen_msg, *ot_sen_msg_1 = this->ot_sen_msg + (this->height-2)*(this->dimension);

        // std::cout<< "   Before opv_gen()" <<std::endl;
        for (int i = 0; i < this->dimension; i++) { // TODO: parallel this 
            opv_gen(this->ggm_master_key[i], this->height, ot_sen_msg_0+i*(this->height-2), ot_sen_msg_1+i*(this->height-2), this->opm_sen_memory+2*i*this->dimension );
        }
        // std::cout<< "   After opv_gen()" <<std::endl;

        // TODO: parallel this
        for (int i = 0; i < this->dimension; i++) {
            for(int j = 0; j < this->dimension; j++) {
                this->column_sum[j] = this->column_sum[j] ^ this->opm_sen_memory[2*i*this->dimension + 2*j];
                this->opm_b[i] = this->opm_b[i] ^ this->opm_sen_memory[2*i*this->dimension + 2*j + 1]; 
                this->opm_a[j] = this->opm_a[j] ^ this->opm_sen_memory[2*i*this->dimension + 2*j + 1];
            }
        }
    }


    void opv_gen(block seed, int height, block *ot_msg_0, block *ot_msg_1, block *opv_memory) {
        //std::cout<< "in opv_gen() 1" <<std::endl;
		TwoKeyPRP *prp = new TwoKeyPRP(zero_block, makeBlock(0, 1));
		prp->node_expand_1to2(opv_memory, seed); //expand `seed` to get two children and store at ggm_tree[0] and ggm_tree[1]
        //std::cout<< "in opv_gen() 2" <<std::endl;
		ot_msg_0[0] = opv_memory[0];
		ot_msg_1[0] = opv_memory[1];
		prp->node_expand_2to4(&opv_memory[0], &opv_memory[0]); /* expand ggm_tree[0] and ggm_tree[1] respectively, then store four extended 
															  children sequentically, starting from ggm_tree[0]. 
															  *ALL CHILDREN NODES ARE STORED IN THE FRONT*.
														    */
        if((height-2) != 1){
            ot_msg_0[1] = opv_memory[0] ^ opv_memory[2];
		    ot_msg_1[1] = opv_memory[1] ^ opv_memory[3];
        }
		for(int h = 2; h < height-1; ++h) {
            if(h != (height-2)) {
			    ot_msg_0[h] = ot_msg_1[h] = zero_block;
            }
			int sz = 1<<h;
            //std::cout<< "\n\nh = " << h << ", sz = " << sz <<std::endl; 
			for(int i = sz-4; i >=0; i-=4) {
				prp->node_expand_4to8(&opv_memory[i*2], &opv_memory[i]);
                if(h != (height-2)) { //no OT messages needed for the last layer, column-wise sum will be sent for the last layer.
                    ot_msg_0[h] ^= opv_memory[i*2];
                    ot_msg_0[h] ^= opv_memory[i*2+2];
                    ot_msg_0[h] ^= opv_memory[i*2+4];
                    ot_msg_0[h] ^= opv_memory[i*2+6];
                    ot_msg_1[h] ^= opv_memory[i*2+1];
                    ot_msg_1[h] ^= opv_memory[i*2+3];
                    ot_msg_1[h] ^= opv_memory[i*2+5];
                    ot_msg_1[h] ^= opv_memory[i*2+7];
                }
                //std::cout<<"opv_memory[i*2+7] = " << i*2+7 << std::endl;
			}
		}
		delete prp;
	}

    void opm_reconstruction() { //`ot_msg_0` and `ot_msg_1` is with length `n*(height-1)`

        bool b[this->height-2];
        for(int i = 0; i < this->dimension; i++) { // TODO: parallel this 
            get_choice_bits(this->permutation[i], b, this->height-2);
            opv_reconstruction(b, this->ot_rev_msg+i*(this->height-2), this->height, opm_rev_memory+i*2*dimension);
            for(int j = 0; j < dimension; j++) {
                // std::cout<<"opm_sen_memory[2*i*dimension+2*j]   = "<< opm_sen_memory[2*i*dimension+2*j] <<  ", opm_rev_memory[2*i*dimension+2*j]   = "<< opm_rev_memory[2*i*dimension+2*j] <<std::endl;
                // std::cout<<"opm_sen_memory[2*i*dimension+2*j+1] = "<< opm_sen_memory[2*i*dimension+2*j+1] <<", opm_rev_memory[2*i*dimension+2*j+1] = "<< opm_rev_memory[2*i*dimension+2*j+1] <<std::endl;
            }
            // std::cout<<std::endl;
        }
        
        // TODO: parallel this
        block *tmp_sum = new block[dimension];
        memset(tmp_sum, 0, sizeof(block)*dimension);

        for (int i = 0; i < dimension; i++) {
            for(int j = 0; j < dimension; j++) {
                tmp_sum[j] = tmp_sum[j] ^ opm_rev_memory[2*i*dimension + 2*j];
>>>>>>> f9de09bd64820fe22cff2ec1e66b1769ab811ef0
            }
        }
        
        // reconstruct the punctured elements 
<<<<<<< HEAD
        for(size_t i = 0; i < dimension; i++) {
            int puncture_idx = permutation[i];
            this->opm_rev_memory[2*i*dimension + 2*puncture_idx] = tmp_sum[puncture_idx] ^ this->column_sum[puncture_idx];  
        }
        delete[] tmp_sum;    
=======
        for(int i = 0; i < dimension; i++) {
            int puncture_idx = permutation[i];
            this->opm_rev_memory[2*i*dimension + 2*puncture_idx] = tmp_sum[puncture_idx] ^ this->column_sum[puncture_idx];  
        }     
        delete[] tmp_sum;   
>>>>>>> f9de09bd64820fe22cff2ec1e66b1769ab811ef0
    }

    void ggm_tree_gen(block seed, int height, block *ot_msg_0, block *ot_msg_1, block *ggm_tree) {
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
<<<<<<< HEAD
		for(size_t h = 2; h < height-1; ++h) {
=======
		for(int h = 2; h < height-1; ++h) {
>>>>>>> f9de09bd64820fe22cff2ec1e66b1769ab811ef0
			ot_msg_0[h] = ot_msg_1[h] = zero_block;
			size_t sz = 1<<h;
			for(int64_t i = sz-4; i >=0; i-=4) {
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

<<<<<<< HEAD
    void ggm_tree_reconstruction(bool *b, block *m, int height, block *ggm_tree) {
		size_t to_fill_idx = 0;
		TwoKeyPRP prp(zero_block, makeBlock(0, 1));
		for(size_t i = 1; i < height; ++i) {
=======

    void ggm_tree_reconstruction(bool *b, block *m, int height, block *ggm_tree) {
		int to_fill_idx = 0;
		TwoKeyPRP prp(zero_block, makeBlock(0, 1));
		for(int i = 1; i < height; ++i) {
>>>>>>> f9de09bd64820fe22cff2ec1e66b1769ab811ef0
			to_fill_idx = to_fill_idx * 2;
			ggm_tree[to_fill_idx] = ggm_tree[to_fill_idx+1] = zero_block; // two undefined children nodes, one will be filled
			if(b[i-1] == false) {// `b[i-1] = 0` <=> `choice_pos[i] = 1`
				layer_recover(height, i, 0, to_fill_idx, m[i-1], &prp, ggm_tree);
				to_fill_idx += 1;
			} else layer_recover(height, i, 1, to_fill_idx+1, m[i-1], &prp, ggm_tree);
		}
	}

<<<<<<< HEAD
	void layer_recover(int height, int depth, size_t lr, size_t to_fill_idx, block sum, TwoKeyPRP *prp, block *ggm_tree) {
=======
	void layer_recover(int height, int depth, int lr, int to_fill_idx, block sum, TwoKeyPRP *prp, block *ggm_tree) {
>>>>>>> f9de09bd64820fe22cff2ec1e66b1769ab811ef0
		int layer_start = 0;
		size_t item_n = 1<<depth;
		block nodes_sum = zero_block;
<<<<<<< HEAD
		size_t lr_start = lr==0?layer_start:(layer_start+1);
		//std::cout<< "\ndepth = " << depth <<std::endl;
		for(size_t i = lr_start; i < item_n; i+=2)
=======
		int lr_start = lr==0?layer_start:(layer_start+1);
		//std::cout<< "\ndepth = " << depth <<std::endl;
		for(int i = lr_start; i < item_n; i+=2)
>>>>>>> f9de09bd64820fe22cff2ec1e66b1769ab811ef0
			nodes_sum = nodes_sum ^ ggm_tree[i];
		ggm_tree[to_fill_idx] = nodes_sum ^ sum;
		if(depth == height-1) return;
		if(item_n == 2)
			prp->node_expand_2to4(&(ggm_tree[0]), &(ggm_tree[0]));
		else {
			for(int i = item_n-4; i >= 0; i-=4){
				prp->node_expand_4to8(&(ggm_tree[i*2]), &(ggm_tree[i]));
            }
		}
	}

    void opv_reconstruction(bool *b, block *m, int height, block *opv_memory) {
		int to_fill_idx = 0;
		TwoKeyPRP prp(zero_block, makeBlock(0, 1));
		for(int i = 1; i < height-1; ++i) { // ends at layer `i = height-2`
			to_fill_idx = to_fill_idx * 2;
			opv_memory[to_fill_idx] = opv_memory[to_fill_idx+1] = zero_block; // two undefined children nodes, one will be filled
			if(b[i-1] == false) {// `b[i-1] = 0` <=> `choice_pos[i] = 1`
				opv_layer_recover(height, i, 0, to_fill_idx, m[i-1], &prp, opv_memory);
				to_fill_idx += 1;
			} else opv_layer_recover(height, i, 1, to_fill_idx+1, m[i-1], &prp, opv_memory);
		}

        // set two zero blocks at layer `height-1`
        to_fill_idx = to_fill_idx * 2;
		opv_memory[to_fill_idx] = opv_memory[to_fill_idx+1] = zero_block;
	}

    // fill layer of `depth` and extend to get layer of `depth+1`
	void opv_layer_recover(int height, int depth, int lr, int to_fill_idx, block sum, TwoKeyPRP *prp, block *ggm_tree) {
		int layer_start = 0;
<<<<<<< HEAD
		size_t item_n = 1<<depth;
		block nodes_sum = zero_block;
		size_t lr_start = lr==0?layer_start:(layer_start+1);
		for(size_t i = lr_start; i < item_n; i+=2)
=======
		int item_n = 1<<depth;
		block nodes_sum = zero_block;
		int lr_start = lr==0?layer_start:(layer_start+1);
		for(int i = lr_start; i < item_n; i+=2)
>>>>>>> f9de09bd64820fe22cff2ec1e66b1769ab811ef0
			nodes_sum = nodes_sum ^ ggm_tree[i];
		ggm_tree[to_fill_idx] = nodes_sum ^ sum;
		//if(depth == height-1) return;
		if(item_n == 2)
			prp->node_expand_2to4(&(ggm_tree[0]), &(ggm_tree[0]));
		else {
<<<<<<< HEAD
			for(size_t i = item_n-4; i >= 0; i-=4){
=======
			for(int i = item_n-4; i >= 0; i-=4){
>>>>>>> f9de09bd64820fe22cff2ec1e66b1769ab811ef0
				prp->node_expand_4to8(&(ggm_tree[i*2]), &(ggm_tree[i]));
            }
		}
	}

<<<<<<< HEAD
    void get_choice_bits(size_t position, bool *b, int length) { // length = bit-lenght(position)
        size_t tmp = position;
		for(size_t i = 0; i < length; i++) {
=======
    void get_choice_bits(int position, bool *b, int length) { // length = bit-lenght(position)
        int tmp = position;
		for(int i = 0; i < length; i++) {
>>>>>>> f9de09bd64820fe22cff2ec1e66b1769ab811ef0
            // std::cout<< "in get_choice_bits for() i = " << i <<std::endl;
			b[length-1-i] = 1 - tmp % 2;
            tmp >>= 1;
		}
	}

<<<<<<< HEAD
    void write_opm_key_to_file(std::string filename, bool is_key = true){
        std::ofstream outfile(filename);
        if(outfile.is_open()) outfile.close();
        else log::error("create a directory to store opm key data\n");
        FileIO fio(filename.c_str(), false);
        fio.send_data(&party, sizeof(int));
        if(party == ALICE) fio.send_data(ggm_master_key, sizeof(block)*dimension);
        else {
            fio.send_data(permutation.data(), sizeof(size_t)*dimension);
            fio.send_data(ot_rev_msg, sizeof(block)*(height-2)*dimension);
        }
    }

    void write_opm_to_file(std::string filename) {
        std::ofstream outfile(filename);
        if(outfile.is_open()) outfile.close();
        else log::error("create a directory to store opm data\n");
        FileIO fio(filename.c_str(), false);
        fio.send_data(&party, sizeof(int));

        if(party == BOB) {
            fio.send_data(permutation.data(), sizeof(size_t)*dimension);
        }

        block *opm_memory = party == ALICE? opm_sen_memory : opm_rev_memory;
        for(size_t i = 0; i < dimension; i++) {
            for(size_t j = 0; j < dimension; j++) {
                fio.send_data(&opm_memory[i*(2*dimension)+2*j+1], sizeof(block));
            }
        }
    }

    void write_opt_to_file(std::string filename) {
        std::ofstream outfile(filename);
        if(outfile.is_open()) outfile.close();
        else log::error("create a directory to store opm data\n");
        FileIO fio(filename.c_str(), false);
        fio.send_data(&party, sizeof(int));

        if(party == BOB) {
            fio.send_data(permutation.data(), sizeof(size_t)*dimension);
            fio.send_data(opm_delta, sizeof(block)*dimension);
        }
        else{
            fio.send_data(this->opm_a, sizeof(block)*dimension);
            fio.send_data(this->opm_b, sizeof(block)*dimension);
        }
    }
    
    void read_opm_key_from_file(std::string filename) {
        FileIO fio(filename.c_str(), true);
        int in_party;
        fio.recv_data(&in_party, sizeof(int));
        if(in_party != party) log::error("wrong party\n");

        if(party == ALICE) {
            fio.recv_data(ggm_master_key, sizeof(block)*dimension);
        }
        else {
            permutation.resize(dimension);
            fio.recv_data(permutation.data(), sizeof(size_t)*dimension);
            fio.recv_data(ot_rev_msg, sizeof(block)*(height-2)*dimension);
        }
    }

    void read_opm_from_file(std::string filename, bool is_key = true){
        FileIO fio(filename.c_str(), true);
        int in_party;
        fio.recv_data(&in_party, sizeof(int));
        if(in_party != party) error("wrong party");

        if(party == BOB) {
            permutation.resize(dimension);
            fio.recv_data(permutation.data(), sizeof(size_t)*dimension);
        }

        block *opm_memory = party == ALICE? opm_sen_memory : opm_rev_memory;
        for(size_t i = 0; i < dimension; i++) {
            for(size_t j = 0; j < dimension; j++) {
                fio.recv_data(&opm_memory[i*(2*dimension)+2*j+1], sizeof(block));
            }
        }
    }

    void read_opt_from_file(std::string filename, bool is_key = true){
        FileIO fio(filename.c_str(), true);
        int in_party;
        fio.recv_data(&in_party, sizeof(int));
        if(in_party != party) log::error("wrong party\n");

        if(party == BOB) {
            permutation.resize(dimension);
            fio.recv_data(permutation.data(), sizeof(size_t)*dimension);
            fio.recv_data(opm_delta, sizeof(block)*dimension);
        }
        else{
            fio.recv_data(this->opm_a, sizeof(block)*dimension);
            fio.recv_data(this->opm_b, sizeof(block)*dimension);
        }
    }

    void test_ggm_tree(int height, size_t index) {
        size_t n = 1<<(height-1);
=======
    template<typename OT>
	void send_opm_f2k(OT *ot, IO *io2, int s) {
		ot->send(ot_sen_msg, ot_sen_msg+height-2, height-2, io2, s);
        io2->send_data(this->column_sum, dimension*sizeof(block));
		io2->send_data(opm_tag, sizeof(block)); // send_data(const void * data, int nbyte). Refer to "emp-tool/io/io_channel.h"
	}

    template<typename OT>
	void recv_opm_f2k(OT *ot, IO *io2, int s) {
		ot->recv(ot_rev_msg, ot_choise_bits, (height-2)*dimension, io2, s);
		io2->recv_data(column_sum, dimension*sizeof(block));
        io2->recv_data(opm_tag, sizeof(block));
	}
    
    
    void test_ggm_tree(int height, int index) {
        int n = 1<<(height-1);
>>>>>>> f9de09bd64820fe22cff2ec1e66b1769ab811ef0
        block seed;
        block *ggm_tree_sen = new block[n];
        block *ggm_tree_rev = new block[n];
        block *ot_msg_sen = new block[2*(height-1)];
        block *ot_msg_rev = new block[height-1];
        bool *b = new bool[height-1];

        prg.random_block(&seed, 1);
        ggm_tree_gen(seed, height, ot_msg_sen, ot_msg_sen+height-1, ggm_tree_sen);
        
        get_choice_bits(index, b, height-1);
<<<<<<< HEAD
        for(size_t i = 0; i < height-1; i++) {// simulating OTs
=======
        for(int i = 0; i < height-1; i++) {// simulating OTs
>>>>>>> f9de09bd64820fe22cff2ec1e66b1769ab811ef0
            ot_msg_rev[i] = (b[i] == false ? ot_msg_sen[i] : ot_msg_sen[height-1+i]);
        }

        ggm_tree_reconstruction(b, ot_msg_rev, height, ggm_tree_rev);

<<<<<<< HEAD
        for(size_t i = 0; i < n; i++) {
=======
        for(int i = 0; i < n; i++) {
>>>>>>> f9de09bd64820fe22cff2ec1e66b1769ab811ef0
            if(i != index && memcmp(ggm_tree_sen+i, ggm_tree_rev+i, sizeof(block)) != 0)
            {
                std::cout<< "ERROR - "<< i << " - " << ggm_tree_sen[i] <<" : "<< ggm_tree_rev[i] <<std::endl;
                //log::error( std::to_string(i) + " " + );
            }
        }

        delete[] ggm_tree_sen;
        delete[] ggm_tree_rev;
        delete[] ot_msg_sen;
        delete[] ot_msg_rev;
        delete[] b;

        std::cout<< "\033[1m\033[32mtest_ggm_tree() SUCCESS" <<std::endl;
    }

    void tes_opv(int height){
        size_t dimension = 1 << (height-2);
        block *opv_sen = new block[2*dimension];
        block *ot_msg_0 = new block[height-2];
        block *ot_msg_1 = new block[height-2];

        block *opv_rev = new block[2*dimension];
        block *ot_rev_msg = new block[height-2];
        bool b[height-2];

        block seed;
        TwoKeyPRP prp(zero_block, makeBlock(0, 1));

        int puncture_indx = 3;

        for(int c = 0; c < 3; c++){
            // simulating opv sender
            prg.random_block(&seed, 1);
            auto t1 = clock_start(); 
            for(size_t i = 0; i < 100; i++) {
                opv_gen(seed, height, ot_msg_0, ot_msg_1, opv_sen);
            }
            std::cout << "\ntime from opv: " << time_from(t1)/100<< " microseconds"  << std::endl;

            // method 2 : use ggm_tree_gen + last-layer extension to generate opv
            // auto t2 = clock_start();
            // for(size_t i = 0; i < 100; i++) {
            //     ggm_tree_gen(seed, height-1, ot_msg_0, ot_msg_1, opv_sen);
            //     for(int i = dimension-1; i >= 0; i--) {
            //         prp.node_expand_1to2(&opv_sen[2*i], opv_sen[i]); 
            //     }
            // }
            // std::cout << "time from ggm_tree+last layer extension: " << time_from(t2)/100<<std::endl;
            
            // simulating opv receiver
            puncture_indx = c;
            get_choice_bits(puncture_indx, b, height-2);
            for(size_t i = 0; i < height-2; i++) {
                ot_rev_msg[i] = b[i] == 0 ? ot_msg_0[i] : ot_msg_1[i];
            }
            opv_reconstruction(b, ot_rev_msg, height, opv_rev);
            for(size_t i = 0; i < 2*dimension; i++) {
                if(i !=2*puncture_indx && i != (2*puncture_indx+1) && !cmpBlock(opv_sen+i, opv_rev+i, 1))
                {
                    std::cout<< "ERROR - "<< i << " - " << opv_sen[i] <<" : "<< opv_rev[i] <<std::endl;
                }
            }
        }
        delete[] opv_sen;
        delete[] ot_msg_0;
        delete[] ot_msg_1;
        delete[] opv_rev;
        delete[] ot_rev_msg;

    }

    void test_OT(int length, int party, int port){ 
	    //NetIO *io = new NetIO(party==ALICE ? nullptr:"127.0.0.1", port);

        // IKNP<NetIO> * iknp = new IKNP<NetIO>(io);
        // std::cout <<"IKNP OT\t"<<double(length)/test_ot<IKNP<NetIO>>(iknp, io, party, length)*1e6<<" OTps"<<std::endl;

        block *msg_0, *msg_1, *rev_msg;;
        msg_0 = new block[length];
        msg_1 = new block[length]; 
        rev_msg = new block[length]; 
        bool *b = new bool[length]; 
        
        this->prg.random_block(msg_0, length);
        this->prg.random_block(msg_1, length);
        this->prg.random_block(rev_msg, length);
        this->prg.random_bool(b, length);

        if(party == ALICE) {
            ot->send(msg_0, msg_1, length); 
        }
        else{
            ot->recv(rev_msg, b, length);
        }
        io->flush();
        std::cout << "test_ot done.\n\n" <<std::endl;

        delete[] msg_0;
        delete[] msg_1;
        delete[] rev_msg;
        delete[] b;
        // delete iknp;
        // delete io;
    }

<<<<<<< HEAD
    void test_opv_with_OT(){
        //int length = (height-2);
        block *ot_sen_msg, *ot_rev_msg, *opv_sen, *opv_rev;
        opv_sen = new block[2*dimension]; 
        opv_rev = new block[2*dimension]; 
        ot_sen_msg = new block[2*(height-2)];
        ot_rev_msg = new block[height-2];
        bool *b = new bool[height-2];

        // FIXME: MUST INITILIZE BLOCKS, OTHERWISE YOU WILL GET ERRORS!
        this->prg.random_block(ot_sen_msg, 2*(height-2));
        this->prg.random_block(ot_rev_msg, (height-2));
        //this->prg.random_bool(b, height-2);

        block seed = zero_block; // only for test
        //prg.random_block(&seed);

        if(party == ALICE){
            this->opv_gen(seed, height, ot_sen_msg, ot_sen_msg+height-2, opv_sen);
            ot->send(ot_sen_msg, ot_sen_msg+height-2, height-2);
        }
        else{
            int puncture_idx = 1;
            get_choice_bits(puncture_idx, b, height-2); 
            ot->recv(ot_rev_msg, b, height-2);
            this->opv_reconstruction(b, ot_rev_msg, height, opv_rev); 
            
            //do correctness check, only for test
            this->opv_gen(seed, height, ot_sen_msg, ot_sen_msg+height-2, opv_sen);
            for(size_t i = 0; i < 2*dimension; i++) {
                if(i !=2*puncture_idx && i != (2*puncture_idx+1) && !cmpBlock(opv_sen+i, opv_rev+i, 1))
                {
                    std::cout<< "ERROR - "<< i << " - " << opv_sen[i] <<" : "<< opv_rev[i] <<std::endl;
                }
            }
        }
        io->flush();
        // std::cout << "test_opv_with_OT() done." <<std::endl;

        delete[] ot_sen_msg;
        delete[] ot_rev_msg;
        delete[] opv_sen; 
        delete[] opv_rev;  
        delete[] b;
    }

    void test_opm() {

        auto t1 = clock_start();
        this->opm_gen();
        std::cout << "time from opm_gen(): " << time_from(t1) << " microseconds" <<std::endl;
       

        // simulating OTs
        block *ot_sen_msg_0 = this->ot_sen_msg, *ot_sen_msg_1 = this->ot_sen_msg + (height-2)*dimension;
        bool b[height-2];
        for(size_t i = 0; i < dimension; i++) {
            get_choice_bits(permutation[i], b, height-2);
            for(size_t j = 0; j < height-2; j++) {
                ot_rev_msg[i*(height-2)+j] = b[j] == false ? ot_sen_msg_0[i*(height-2)+j] : ot_sen_msg_1[i*(height-2)+j]; // NOTE: b[j] not b[i] !!!
            }
        }

        auto t2 = clock_start(); 
        this->opm_reconstruction();
        std::cout << "time from opm_reconstruction(): " << time_from(t2) << " microseconds" <<std::endl;

        //simulating consistency check
        block *chi = new block[dimension*dimension]; 
		Hash hash;
		block digest[2];

        auto t3 = clock_start();
        // consistency check method 1
		hash.hash_once(digest, this->column_sum, dimension*sizeof(block));
		uni_hash_coeff_gen(chi, digest[0], dimension*dimension);
        for(size_t i = 0; i < dimension; i++) {
            std::cout << " column sum i = " << i << " : " << column_sum[i];
        }
        std::cout<<std::endl; 
        std::cout << "digest[0]: " << digest[0] <<std::endl;
        std::cout << "chi[0]: " << chi[0] <<std::endl;
        
        block r_sen = zero_block, r_rev = zero_block;
		block r1, r2;
        for(size_t i = 0; i < dimension; i++) {
            for(size_t j = 0; j < dimension; j++) {
                gfmul(this->opm_sen_memory[2*i*dimension+2*j], chi[i*dimension+j], &r1);
                gfmul(this->opm_rev_memory[2*i*dimension+2*j], chi[i*dimension+j], &r2);
                r_sen = r_sen ^ r1;
                r_rev = r_rev ^ r2;
            }
        }
        if(!cmpBlock(&r_sen, &r_rev, 1)) 
        {
            log::error("OPM check fails!\n");
        }
        std::cout << "r_sen: " << r_sen << ", r_rev: " << r_rev <<std::endl;
        std::cout << "time from OPM check method 1 : " << time_from(t3) << " microseconds" <<std::endl;

        // consistency check method 2
        // auto t4 = clock_start(); 
        // block digest_rev[2];
        // hash.hash_once(digest, opm_sen_memory, 2*dimension*dimension*sizeof(block));
        // hash.hash_once(digest_rev, opm_sen_memory, 2*dimension*dimension*sizeof(block));
        // if(!cmpBlock(digest, digest_rev, 2)) {
        //     std::cout<<"OPM check fails."<<std::endl;
        // }
        // std::cout << "time from OPM check method 2: " << time_from(t4) << " microseconds" <<std::endl;
		
        // check OPT relation
        block *opt_check = new block[dimension];
        for(size_t i = 0; i < dimension; i++) {
            // delta = p(a)+b 
            opt_check[i] = opm_delta[i] ^ opm_a[permutation[i]] ^ opm_b[i];
            if(!cmpBlock(opt_check+i, &zero_block, 1)) {
                log::error("opt check fails.\n");
            }
        }
        delete[] opt_check;  

        delete[] chi;
    }

    void test_opm_with_OT() {
        block *chi = new block[dimension*dimension]; 
        Hash hash;
        block digest[2];
        block sen_tag = zero_block, rev_tag = zero_block;

        if(party == ALICE){
            //auto t1 = clock_start(); 
            this->opm_gen();
            //std::cout << "time from opm_gen(): " << time_from(t1) << " microseconds" <<std::endl;

            hash.hash_once(digest, column_sum, dimension*sizeof(block)); //`column_sum` was updated by function `opm_gen()` 
            uni_hash_coeff_gen(chi, digest[0], dimension*dimension);
            block r1;
            for(size_t i = 0; i < dimension; i++) {
                for(size_t j = 0; j < dimension; j++) {
                    gfmul(this->opm_sen_memory[2*i*dimension+2*j], chi[i*dimension+j], &r1);
                    sen_tag = sen_tag ^ r1;
                }
            }

            ot->send(ot_sen_msg, ot_sen_msg+(height-2)*dimension, (height-2)*dimension);
            //std::cout<< "ot sen msg: " << ot_sen_msg[0] << " " << ot_sen_msg[(height-2)*dimension] <<std::endl;
            io->send_data(column_sum, dimension*sizeof(block));io->flush();
            //std::cout<<"sen column_sum: "<< *column_sum <<std::endl;
            hash.hash_once(digest, &sen_tag, sizeof(block));
            io->send_data(digest, 2*sizeof(block));io->flush();
            //std::cout<<"sen digest: "<< digest[0] << " " << digest[1] <<std::endl;
            io->recv_data(&rev_tag, sizeof(block));//io->flush();
            // do check
            if(!cmpBlock(&rev_tag, &sen_tag, 1)) 
            {
                log::error("hash from receiver does not match sender's tag\n");
            }
        }
        else{
            for(size_t i = 0; i < dimension; i++) {
                get_choice_bits(permutation[i], ot_choise_bits+i*(height-2), height-2);
            }

            ot->recv(ot_rev_msg, ot_choise_bits, (height-2)*dimension);
            //std::cout<< "ot rev msg: " << *ot_rev_msg <<std::endl;
            io->recv_data(column_sum, dimension*sizeof(block));
            
            //io->flush();

            auto t2 = clock_start(); 
            this->opm_reconstruction();
            //std::cout << "time from opm_reconstruction(): " << time_from(t2) << " microseconds" <<std::endl;
            //std::cout<<"rev column_sum: "<< *column_sum <<std::endl;
            // consistency check
            Hash hash;
            block digest[2];

            auto t3 = clock_start();
            hash.hash_once(digest, column_sum, dimension*sizeof(block));
            uni_hash_coeff_gen(chi, digest[0], dimension*dimension);
            block r1;
            for(size_t i = 0; i < dimension; i++) {
                for(size_t j = 0; j < dimension; j++) {
                    gfmul(this->opm_rev_memory[2*i*dimension+2*j], chi[i*dimension+j], &r1);
                    rev_tag = rev_tag ^ r1;
                }
            }

            block hash_sen[2];
            io->recv_data(hash_sen, 2*sizeof(block));
            hash.hash_once(digest, &rev_tag, sizeof(block));
            //std::cout<<"rev digest: "<< digest[0] << " " << digest[1] <<std::endl;
            if(!cmpBlock(hash_sen, digest, 2)) 
            {
                log::error("hash from sender does not match receiver's tag\n");
                //return ;//FIXME: 
            }

            // send back the tag for check
            io->send_data(&rev_tag, sizeof(block));io->flush();

            //std::cout << "time from OPM check method 1 : " << time_from(t3) << " microseconds" <<std::endl;
        }
        //io->flush();

        delete[] chi;
        //std::cout << "test_opm_with_OT DONE." <<std::endl;
=======
        delete[] ggm_tree_sen;
        delete[] ggm_tree_rev;
        delete[] ot_msg_sen;
        delete[] ot_msg_rev;
        delete[] b;

        std::cout<< "\033[1m\033[32mtest_ggm_tree() SUCCESS" <<std::endl;
    }


    void tes_opv(int height){
        int dimension = 1 << (height-2);
        block *opv_sen = new block[2*dimension];
        block *ot_msg_0 = new block[height-2];
        block *ot_msg_1 = new block[height-2];

        block *opv_rev = new block[2*dimension];
        block *ot_rev_msg = new block[height-2];
        bool b[height-2];

        block seed;
        TwoKeyPRP prp(zero_block, makeBlock(0, 1));

        int puncture_indx = 3;

        for(int c = 0; c < 3; c++){
            // simulating opv sender
            prg.random_block(&seed, 1);
            auto t1 = clock_start(); 
            for(int i = 0; i < 100; i++) {
                opv_gen(seed, height, ot_msg_0, ot_msg_1, opv_sen);
            }
            std::cout << "\ntime from opv: " << time_from(t1)/100<< " microseconds"  << std::endl;

            // method 2 : use ggm_tree_gen + last-layer extension to generate opv
            // auto t2 = clock_start();
            // for(int i = 0; i < 100; i++) {
            //     ggm_tree_gen(seed, height-1, ot_msg_0, ot_msg_1, opv_sen);
            //     for(int i = dimension-1; i >= 0; i--) {
            //         prp.node_expand_1to2(&opv_sen[2*i], opv_sen[i]); 
            //     }
            // }
            // std::cout << "time from ggm_tree+last layer extension: " << time_from(t2)/100<<std::endl;
            
            // simulating opv receiver
            puncture_indx = c;
            get_choice_bits(puncture_indx, b, height-2);
            for(int i = 0; i < height-2; i++) {
                ot_rev_msg[i] = b[i] == 0 ? ot_msg_0[i] : ot_msg_1[i];
            }
            opv_reconstruction(b, ot_rev_msg, height, opv_rev);
            for(int i = 0; i < 2*dimension; i++) {
                if(i !=2*puncture_indx && i != (2*puncture_indx+1) && !cmpBlock(opv_sen+i, opv_rev+i, 1))
                {
                    std::cout<< "ERROR - "<< i << " - " << opv_sen[i] <<" : "<< opv_rev[i] <<std::endl;
                }
            }
        }
        delete[] opv_sen;
        delete[] ot_msg_0;
        delete[] ot_msg_1;
        delete[] opv_rev;
        delete[] ot_rev_msg;

    }

    void test_opm() {

        auto t1 = clock_start(); 
        this->opm_gen();
        std::cout << "time from opm_gen(): " << time_from(t1) << " microseconds" <<std::endl;
       

        // simulating OTs
        block *ot_sen_msg_0 = this->ot_sen_msg, *ot_sen_msg_1 = this->ot_sen_msg + (height-2)*dimension;
        bool b[height-2];
        for(int i = 0; i < dimension; i++) {
            get_choice_bits(permutation[i], b, height-2);
            for(int j = 0; j < height-2; j++) {
                ot_rev_msg[i*(height-2)+j] = b[j] == false ? ot_sen_msg_0[i*(height-2)+j] : ot_sen_msg_1[i*(height-2)+j]; // NOTE: b[j] not b[i] !!!
            }
        }

        auto t2 = clock_start(); 
        this->opm_reconstruction();
        std::cout << "time from opm_reconstruction(): " << time_from(t2) << " microseconds" <<std::endl;

        //simulating consistency check
        block *chi = new block[dimension*dimension]; 
		Hash hash;
		block digest[2];

        auto t3 = clock_start();
        // consistency check method 1
		hash.hash_once(digest, this->column_sum, dimension*sizeof(block));
		uni_hash_coeff_gen(chi, digest[0], dimension*dimension);
        
        block r_sen = zero_block, r_rev = zero_block;
		block r1, r2;
        for(int i = 0; i < dimension; i++) {
            for(int j = 0; j < dimension; j++) {
                gfmul(this->opm_sen_memory[2*i*dimension+2*j], chi[i*dimension+j], &r1);
                gfmul(this->opm_rev_memory[2*i*dimension+2*j], chi[i*dimension+j], &r2);
                r_sen = r_sen ^ r1;
                r_rev = r_rev ^ r2;
            }
        }
        if(!cmpBlock(&r_sen, &r_rev, 1)) 
        {
            std::cout<<"OPM check fails."<<std::endl;
        }
        std::cout << "time from OPM check method 1 : " << time_from(t3) << " microseconds" <<std::endl;

        // consistency check method 2
        // auto t4 = clock_start(); 
        // block digest_rev[2];
        // hash.hash_once(digest, opm_sen_memory, 2*dimension*dimension*sizeof(block));
        // hash.hash_once(digest_rev, opm_sen_memory, 2*dimension*dimension*sizeof(block));
        // if(!cmpBlock(digest, digest_rev, 2)) {
        //     std::cout<<"OPM check fails."<<std::endl;
        // }
        // std::cout << "time from OPM check method 2: " << time_from(t4) << " microseconds" <<std::endl;
		
        delete[] chi;
>>>>>>> f9de09bd64820fe22cff2ec1e66b1769ab811ef0
    }
};

#endif 