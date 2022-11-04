#ifndef SPCOT_SENDER_H__
#define SPCOT_SENDER_H__
#include <iostream>
#include "emp-tool/emp-tool.h"
#include "emp-ot/emp-ot.h"
#include "emp-ot/ferret/twokeyprp.h"

using namespace emp;

template<typename IO>
class SPCOT_Sender { public:
	block seed;
	block delta;
	block *ggm_tree, *m;
	IO *io;
	int depth, leave_n; // `leave_n` denotes number of leaves of the GGM tree, ie, number of nodes of layer `depth-1` 
	PRG prg;
	block secret_sum_f2;

	SPCOT_Sender(IO *io, int depth_in) {
		initialization(io, depth_in);
		prg.random_block(&seed, 1);
	}

	void initialization(IO *io, int depth_in) {
		this->io = io;
		this->depth = depth_in;
		this->leave_n = 1<<(this->depth-1); // depth: the number of nodes from the root to a leaf. 
		m = new block[(depth-1)*2]; /* `m` is used for sender's OT messsages computed from GGM tree. 
								 		Each level produces two messages from level 2 to depth, thus 
										needing (depth-1)*2 in total
									*/
	}

	~SPCOT_Sender() {
		delete[] m;
	}

	// generate GGM tree, transfer secret, F2^k
	void compute(block* ggm_tree_mem, block secret) {
		this->delta = secret;
		ggm_tree_gen(m, m+depth-1, ggm_tree_mem, secret); // `m+depth-1` is the entry address for the second half of OT messages  
	}

	// send the nodes by oblivious transfer, F2^k
	template<typename OT>
	void send_f2k(OT * ot, IO * io2, int s) {
		ot->send(m, &m[depth-1], depth-1, io2, s); /* void send(const block* data0, const block* data1, int64_t length): 
														`data0`	: the first half of OT messages;
														`data1`	: the second half of OT messages;
														`length`: number of pairs to send.									
													*/
		io2->send_data(&secret_sum_f2, sizeof(block)); // send_data(const void * data, int nbyte). Refer to "emp-tool/io/io_channel.h"
	}

	void ggm_tree_gen(block *ot_msg_0, block *ot_msg_1, block* ggm_tree_mem, block secret) { 
		ggm_tree_gen(ot_msg_0, ot_msg_1, ggm_tree_mem); // `ggm_tree_mem` is passed from some outter caller
		secret_sum_f2 = zero_block;
		block one = makeBlock(0xFFFFFFFFFFFFFFFFLL,0xFFFFFFFFFFFFFFFELL); 
		for(int i = 0; i < leave_n; ++i) { // all leaves are stored in the front of ggm_tree, see function ggm_tree_gen for details. 
			ggm_tree[i] = ggm_tree[i] & one; // TODO: why & ? 
			secret_sum_f2 = secret_sum_f2 ^ ggm_tree[i];
		}
		secret_sum_f2 = secret_sum_f2 ^ secret;
	}

	// generate GGM tree from the top
	void ggm_tree_gen(block *ot_msg_0, block *ot_msg_1, block* ggm_tree_mem) {
		this->ggm_tree = ggm_tree_mem; // this->ggm_tree points to `ggm_tree_mem` from the outter caller
		TwoKeyPRP *prp = new TwoKeyPRP(zero_block, makeBlock(0, 1));
		prp->node_expand_1to2(ggm_tree, seed); // expand `seed` to get two children and store at ggm_tree[0] and ggm_tree[1]
		ot_msg_0[0] = ggm_tree[0];
		ot_msg_1[0] = ggm_tree[1];
		prp->node_expand_2to4(&ggm_tree[0], &ggm_tree[0]); /* expand ggm_tree[0] and ggm_tree[1] respectively, then store four extended 
															  children sequentically, beginning from ggm_tree[0]. 
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

	void consistency_check_msg_gen(block *V) {
		// X
		block *chi = new block[leave_n];
		Hash hash;
		block digest[2];
		hash.hash_once(digest, &secret_sum_f2, sizeof(block)); // `digest = hash(secret_sum_f2)`
		uni_hash_coeff_gen(chi, digest[0], leave_n); // use `digest[0]` to generate `leave_n` coefficients and store at `chi`, which are [digest[0], digest[0]^2, ..., digest[0]^{leave_n-1}]

		vector_inn_prdt_sum_red(V, chi, ggm_tree, leave_n); // `V = Sum_i chi[i] * ggm_tree[i]`
		delete[] chi;
	}
};

#endif
