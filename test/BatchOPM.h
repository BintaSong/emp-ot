#ifndef BATCH_OPM_H__
#define BATCH_OPM_H__

#include <emp-tool/emp-tool.h>
#include "emp-ot/emp-ot.h"
#include "test/shuffle_common.h"
#include <iostream>
#include <algorithm>
#include <random>

#define MALICIOUS

using namespace emp;
using std::future;

template <typename IO>
class BatchOPM
{
public:
    int party, n_threads;
    size_t depth, height, dimension;

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
    bool *ot_choise_bits = nullptr;

    ThreadPool *pool;

    PRG prg;

    BatchOPM(int party, IO *io, OT<IO> *ot, int threads, size_t depth, bool pre_setup = false)
    { // height > 3 for OPM
        this->party = party;
        this->io = io;
        this->ot = ot;
        this->n_threads = threads;
        this->depth = depth;
        this->height = depth + 1;
#ifdef MALICIOUS
        this->height = depth + 2;
#endif
        this->dimension = 1 << (depth);
        this->column_sum = new block[dimension];
        memset(this->column_sum, 0, this->dimension * sizeof(block));
        this->opm_tag = new block;
        // if (party == BOB) receiver_init();
        // else sender_init();

        pool = new ThreadPool(threads);

        receiver_init(pre_setup);
        sender_init(pre_setup);
    }

    void receiver_init(bool pre_setup = false)
    {
        // std::cout<< "new opm_rev_memory" <<std::endl;
        this->opm_rev_memory = new block[2 * this->dimension * this->dimension];
        memset(this->opm_rev_memory, 0, 2 * this->dimension * this->dimension * sizeof(block));

        this->opm_delta = new block[this->dimension];
        memset(this->opm_delta, 0, this->dimension * sizeof(block));

        ot_rev_msg = new block[(height - 2) * dimension];
        ot_choise_bits = new bool[(height - 2) * dimension];

        if (pre_setup)
        { // setup receiver from local stored file

            return;
        }

        this->random_permutation();
    }

    void sender_init(bool pre_setup = false)
    {

        this->ggm_master_key = new block[this->dimension];
        prg.random_block(this->ggm_master_key, this->dimension);

        this->opm_sen_memory = new block[2 * this->dimension * this->dimension];
        this->opm_a = new block[this->dimension];
        this->opm_b = new block[this->dimension];
        memset(this->opm_sen_memory, 0, 2 * this->dimension * this->dimension * sizeof(block));
        memset(this->opm_a, 0, this->dimension * sizeof(block));
        memset(this->opm_b, 0, this->dimension * sizeof(block));
        this->ot_sen_msg = new block[2 * (this->height - 2) * this->dimension]; // layers [2, height-1] need OT messages
        if (pre_setup)
        { // setup sender from local stored file

            return;
        }
    }

    ~BatchOPM()
    {
        if (ggm_master_key != nullptr)
            delete[] ggm_master_key;

        delete[] column_sum;
        delete opm_tag;

        if (opm_sen_memory != nullptr)
            delete[] opm_sen_memory;
        if (opm_rev_memory != nullptr)
            delete[] opm_rev_memory;

        if (opm_a != nullptr)
            delete[] opm_a;
        if (opm_b != nullptr)
            delete[] opm_b;
        if (opm_delta != nullptr)
            delete[] opm_delta;
        if (ot_sen_msg != nullptr)
            delete[] ot_sen_msg;
        if (ot_rev_msg != nullptr)
            delete[] ot_rev_msg;
        if (ot_choise_bits != nullptr)
            delete[] ot_choise_bits;
        delete pool;
        // std::cout<<"~OPM() end."<<std::endl;
    }

    // TODO: use cryptigraphic secure random source
    void random_permutation()
    {
        for (int i = 0; i < this->dimension; i++)
        {
            this->permutation.push_back(i);
        }

        std::random_device rand_div("/dev/urandom");
        auto rng = std::default_random_engine{rand_div()};
        std::shuffle(permutation.begin(), permutation.end(), rng);

        permutation_inverse.resize(dimension);
        for (size_t i = 0; i < dimension; i++)
        {
            permutation_inverse[permutation[i]] = i;
        }
    }

    void get_random_permutation(const size_t dimension, std::vector<size_t> &permutation, std::vector<size_t> &permutation_inverse)
    {
        permutation.resize(0);
        for (size_t i = 0; i < dimension; i++)
        {
            permutation.push_back(i);
        }

        std::random_device rand_div("/dev/urandom");
        auto rng = std::default_random_engine{rand_div()};
        std::shuffle(permutation.begin(), permutation.end(), rng);

        permutation_inverse.resize(dimension);
        for (size_t i = 0; i < dimension; i++)
        {
            permutation_inverse[permutation[i]] = i;
        }
    }

    void opm_gen(const int height, const size_t dimension, const block *ggm_master_key, block *opm_sen_memory, block *ot_sen_msg_0, block *ot_sen_msg_1, block *opm_a, block *opm_b, block *column_sum)
    { //`ot_msg_0` and `ot_msg_1` is with length `n*(height-1)`

        // block *ot_sen_msg_0 = ot_sen_msg, *ot_sen_msg_1 = ot_sen_msg + (height-2)*dimension;

        // std::cout<< "   Before opv_gen()" <<std::endl;
        for (size_t i = 0; i < dimension; i++)
        { // TODO: parallel this
            opv_gen(ggm_master_key[i], height, ot_sen_msg_0 + i * (height - 2), ot_sen_msg_1 + i * (height - 2), opm_sen_memory + 2 * i * dimension);
        }
        // std::cout<< "   After opv_gen()" <<std::endl;

        // TODO: parallel this
        memset(column_sum, 0, dimension * sizeof(block));

        for (size_t i = 0; i < dimension; i++)
        {
            for (int j = 0; j < dimension; j++)
            {
                column_sum[j] = column_sum[j] ^ opm_sen_memory[2 * i * dimension + 2 * j];
                opm_b[i] = opm_b[i] ^ opm_sen_memory[2 * i * dimension + 2 * j + 1];
                opm_a[j] = opm_a[j] ^ opm_sen_memory[2 * i * dimension + 2 * j + 1];
            }
        }
    }

    void opv_gen(const block seed, const size_t height, block *ot_msg_0, block *ot_msg_1, block *opv_memory)
    {
        // std::cout<< "in opv_gen() 1" <<std::endl;
        TwoKeyPRP *prp = new TwoKeyPRP(zero_block, makeBlock(0, 1));
        prp->node_expand_1to2(opv_memory, seed); // expand `seed` to get two children and store at ggm_tree[0] and ggm_tree[1]
        // std::cout<< "in opv_gen() 2" <<std::endl;
        ot_msg_0[0] = opv_memory[0];
        ot_msg_1[0] = opv_memory[1];
        prp->node_expand_2to4(&opv_memory[0], &opv_memory[0]); /* expand ggm_tree[0] and ggm_tree[1] respectively, then store four extended
                                                              children sequentically, starting from ggm_tree[0].
                                                              *ALL CHILDREN NODES ARE STORED IN THE FRONT*.
                                                            */
        if ((height - 2) != 1)
        {
            ot_msg_0[1] = opv_memory[0] ^ opv_memory[2];
            ot_msg_1[1] = opv_memory[1] ^ opv_memory[3];
        }
        for (int h = 2; h < height - 1; ++h)
        {
            if (h != (height - 2))
            {
                ot_msg_0[h] = ot_msg_1[h] = zero_block;
            }
            int sz = 1 << h;
            // std::cout<< "\n\nh = " << h << ", sz = " << sz <<std::endl;
            for (int i = sz - 4; i >= 0; i -= 4)
            {
                prp->node_expand_4to8(&opv_memory[i * 2], &opv_memory[i]);
                if (h != (height - 2))
                { // no OT messages needed for the last layer, column-wise sum will be sent for the last layer.
                    ot_msg_0[h] ^= opv_memory[i * 2];
                    ot_msg_0[h] ^= opv_memory[i * 2 + 2];
                    ot_msg_0[h] ^= opv_memory[i * 2 + 4];
                    ot_msg_0[h] ^= opv_memory[i * 2 + 6];
                    ot_msg_1[h] ^= opv_memory[i * 2 + 1];
                    ot_msg_1[h] ^= opv_memory[i * 2 + 3];
                    ot_msg_1[h] ^= opv_memory[i * 2 + 5];
                    ot_msg_1[h] ^= opv_memory[i * 2 + 7];
                }
                // std::cout<<"opv_memory[i*2+7] = " << i*2+7 << std::endl;
            }
        }
        delete prp;
    }

    void opm_reconstruction(const int height, const size_t dimension, const std::vector<size_t> &permutation, const std::vector<size_t> permutation_inverse, const block *ot_rev_msg, const block *column_sum, block *opm_rev_memory, block *opm_delta)
    { //`ot_msg_0` and `ot_msg_1` is with length `n*(height-1)`

        bool b[height - 2];
        for (size_t i = 0; i < dimension; i++)
        { // TODO: parallel this
            get_choice_bits(permutation[i], b, this->height - 2);
            opv_reconstruction(b, ot_rev_msg + i * (height - 2), height, opm_rev_memory + i * 2 * dimension);
            // log::green("\n\nnew recon:\n");
            // for(int j = 0; j < dimension; j++) {
            //     std::cout<<"opm_sen_memory[2*i*dimension+2*j]   = "<< opm_sen_memory[2*i*dimension+2*j] <<  ", opm_rev_memory[2*i*dimension+2*j]   = "<< opm_rev_memory[2*i*dimension+2*j] <<std::endl;
            //     std::cout<<"opm_sen_memory[2*i*dimension+2*j+1] = "<< opm_sen_memory[2*i*dimension+2*j+1] <<", opm_rev_memory[2*i*dimension+2*j+1] = "<< opm_rev_memory[2*i*dimension+2*j+1] <<std::endl;
            // }
            // std::cout<<std::endl;
        }

        // update opm_delta and prepare for opm sacrifice TODO: parallel this

        block *tmp_sum = new block[dimension];
        memset(tmp_sum, 0, sizeof(block) * dimension);

        for (size_t i = 0; i < dimension; i++)
        {
            for (int j = 0; j < dimension; j++)
            {

                tmp_sum[j] = tmp_sum[j] ^ opm_rev_memory[2 * i * dimension + 2 * j];

                if (j != permutation[i])
                {
                    opm_delta[i] ^= opm_rev_memory[2 * i * dimension + 2 * j + 1];
                    opm_delta[permutation_inverse[j]] ^= opm_rev_memory[2 * i * dimension + 2 * j + 1];
                }
            }
        }

        // reconstruct the punctured elements
        for (size_t i = 0; i < dimension; i++)
        {
            int puncture_idx = permutation[i];
            opm_rev_memory[2 * i * dimension + 2 * puncture_idx] = tmp_sum[puncture_idx] ^ column_sum[puncture_idx];
        }
        delete[] tmp_sum;
    }

    void opv_reconstruction(const bool *b, const block *m, const int height, block *opv_memory)
    { //`m` is the received ot msgs
        int to_fill_idx = 0;
        TwoKeyPRP prp(zero_block, makeBlock(0, 1));
        for (int i = 1; i < height - 1; ++i)
        { // ends at layer `i = height-2`
            to_fill_idx = to_fill_idx * 2;
            opv_memory[to_fill_idx] = opv_memory[to_fill_idx + 1] = zero_block; // two undefined children nodes, one will be filled
            if (b[i - 1] == false)
            { // `b[i-1] = 0` <=> `choice_pos[i] = 1`
                opv_layer_recover(height, i, 0, to_fill_idx, m[i - 1], &prp, opv_memory);
                to_fill_idx += 1;
            }
            else
                opv_layer_recover(height, i, 1, to_fill_idx + 1, m[i - 1], &prp, opv_memory);
        }

        // set two zero blocks at layer `height-1`
        to_fill_idx = to_fill_idx * 2;
        opv_memory[to_fill_idx] = opv_memory[to_fill_idx + 1] = zero_block;
    }

    void ggm_tree_gen(block seed, int height, block *ot_msg_0, block *ot_msg_1, block *ggm_tree)
    {
        // this->ggm_tree = ggm_tree_mem; // this->ggm_tree points to `ggm_tree_mem` from the outter caller
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
        for (size_t h = 2; h < height - 1; ++h)
        {
            ot_msg_0[h] = ot_msg_1[h] = zero_block;
            size_t sz = 1 << h;
            for (int64_t i = sz - 4; i >= 0; i -= 4) // FIXME: DO NOT USE UNSIGNED INT !
            {
                prp->node_expand_4to8(&ggm_tree[i * 2], &ggm_tree[i]);
                ot_msg_0[h] ^= ggm_tree[i * 2];
                ot_msg_0[h] ^= ggm_tree[i * 2 + 2];
                ot_msg_0[h] ^= ggm_tree[i * 2 + 4];
                ot_msg_0[h] ^= ggm_tree[i * 2 + 6];
                ot_msg_1[h] ^= ggm_tree[i * 2 + 1];
                ot_msg_1[h] ^= ggm_tree[i * 2 + 3];
                ot_msg_1[h] ^= ggm_tree[i * 2 + 5];
                ot_msg_1[h] ^= ggm_tree[i * 2 + 7];
            }
        }
        delete prp;
    }

    void ggm_tree_reconstruction(bool *b, block *m, int height, block *ggm_tree)
    {
        size_t to_fill_idx = 0;
        TwoKeyPRP prp(zero_block, makeBlock(0, 1));
        for (size_t i = 1; i < height; ++i)
        {
            to_fill_idx = to_fill_idx * 2;
            ggm_tree[to_fill_idx] = ggm_tree[to_fill_idx + 1] = zero_block; // two undefined children nodes, one will be filled
            if (b[i - 1] == false)
            { // `b[i-1] = 0` <=> `choice_pos[i] = 1`
                layer_recover(height, i, 0, to_fill_idx, m[i - 1], &prp, ggm_tree);
                to_fill_idx += 1;
            }
            else
                layer_recover(height, i, 1, to_fill_idx + 1, m[i - 1], &prp, ggm_tree);
        }
    }

    void layer_recover(int height, int depth, size_t lr, size_t to_fill_idx, block sum, TwoKeyPRP *prp, block *ggm_tree)
    {
        int layer_start = 0;
        size_t item_n = 1 << depth;
        block nodes_sum = zero_block;
        size_t lr_start = lr == 0 ? layer_start : (layer_start + 1);
        // std::cout<< "\ndepth = " << depth <<std::endl;
        for (size_t i = lr_start; i < item_n; i += 2)
            nodes_sum = nodes_sum ^ ggm_tree[i];
        ggm_tree[to_fill_idx] = nodes_sum ^ sum;
        if (depth == height - 1)
            return;
        if (item_n == 2)
            prp->node_expand_2to4(&(ggm_tree[0]), &(ggm_tree[0]));
        else
        {
            for (int i = item_n - 4; i >= 0; i -= 4)
            {
                prp->node_expand_4to8(&(ggm_tree[i * 2]), &(ggm_tree[i]));
            }
        }
    }

    // fill layer of `depth` and extend to get layer of `depth+1`
    void opv_layer_recover(int height, int depth, int lr, int to_fill_idx, block sum, TwoKeyPRP *prp, block *ggm_tree)
    {
        int layer_start = 0;
        int item_n = 1 << depth;
        block nodes_sum = zero_block;
        int lr_start = lr == 0 ? layer_start : (layer_start + 1);
        for (int i = lr_start; i < item_n; i += 2)
            nodes_sum = nodes_sum ^ ggm_tree[i];
        ggm_tree[to_fill_idx] = nodes_sum ^ sum;
        // if(depth == height-1) return;
        if (item_n == 2)
            prp->node_expand_2to4(&(ggm_tree[0]), &(ggm_tree[0]));
        else
        {
            for (int i = item_n - 4; i >= 0; i -= 4)
            {
                prp->node_expand_4to8(&(ggm_tree[i * 2]), &(ggm_tree[i]));
            }
        }
    }

    void get_choice_bits(size_t position, bool *b, size_t length)
    { // length = bit-lenght(position)
        size_t tmp = position;
        for (int i = 0; i < length; i++)
        {
            // std::cout<< "in get_choice_bits for() i = " << i <<std::endl;
            b[length - 1 - i] = 1 - tmp % 2;
            tmp >>= 1;
        }
    }

    void write_opm_key_to_file(std::string filename, bool is_key = true)
    {
        std::ofstream outfile(filename);
        if (outfile.is_open())
            outfile.close();
        else
            log::error("create a directory to store opm key data\n");
        FileIO fio(filename.c_str(), false);
        fio.send_data(&party, sizeof(int));
        if (party == ALICE)
            fio.send_data(ggm_master_key, sizeof(block) * dimension);
        else
        {
            fio.send_data(permutation.data(), sizeof(int) * dimension);
            fio.send_data(ot_rev_msg, sizeof(block) * (height - 2) * dimension);
        }
    }

    void write_opm_to_file(std::string filename)
    {
        std::ofstream outfile(filename);
        if (outfile.is_open())
            outfile.close();
        else
            log::error("create a directory to store opm data\n");
        FileIO fio(filename.c_str(), false);
        fio.send_data(&party, sizeof(int));

        if (party == BOB)
        {
            fio.send_data(permutation.data(), sizeof(int) * dimension);
        }

        block *opm_memory = party == ALICE ? opm_sen_memory : opm_rev_memory;
        for (size_t i = 0; i < dimension; i++)
        {
            for (int j = 0; j < dimension; j++)
            {
                fio.send_data(&opm_memory[i * (2 * dimension) + 2 * j + 1], sizeof(block));
            }
        }
    }

    void read_opm_key_from_file(std::string filename)
    {
        FileIO fio(filename.c_str(), true);
        int in_party;
        fio.recv_data(&in_party, sizeof(int));
        if (in_party != party)
            log::error("wrong party\n");

        if (party == ALICE)
        {
            fio.recv_data(ggm_master_key, sizeof(block) * dimension);
        }
        else
        {
            permutation.resize(dimension);
            fio.recv_data(permutation.data(), sizeof(int) * dimension);
            fio.recv_data(ot_rev_msg, sizeof(block) * (height - 2) * dimension);
        }
    }

    void read_opm_from_file(std::string filename, bool is_key = true)
    {
        FileIO fio(filename.c_str(), true);
        int in_party;
        fio.recv_data(&in_party, sizeof(int));
        if (in_party != party)
            error("wrong party");

        if (party == BOB)
        {
            permutation.resize(dimension);
            fio.recv_data(permutation.data(), sizeof(int) * dimension);
        }

        block *opm_memory = party == ALICE ? opm_sen_memory : opm_rev_memory;
        for (size_t i = 0; i < dimension; i++)
        {
            for (int j = 0; j < dimension; j++)
            {
                fio.recv_data(&opm_memory[i * (2 * dimension) + 2 * j + 1], sizeof(block));
            }
        }
    }

    void read_opt_from_file(std::string filename, bool is_key = true)
    {
        FileIO fio(filename.c_str(), true);
        int in_party;
        fio.recv_data(&in_party, sizeof(int));
        if (in_party != party)
            log::error("wrong party\n");

        if (party == BOB)
        {
            permutation.resize(dimension);
            fio.recv_data(permutation.data(), sizeof(int) * dimension);
            fio.recv_data(opm_delta, sizeof(block) * dimension);
        }
        else
        {
            fio.recv_data(this->opm_a, sizeof(block) * dimension);
            fio.recv_data(this->opm_b, sizeof(block) * dimension);
        }
    }

    void tes_opv(int height)
    {
        size_t dimension = 1 << (height - 2);
        block *opv_sen = new block[2 * dimension];
        block *ot_msg_0 = new block[height - 2];
        block *ot_msg_1 = new block[height - 2];

        block *opv_rev = new block[2 * dimension];
        block *ot_rev_msg = new block[height - 2];
        bool b[height - 2];

        block seed;
        TwoKeyPRP prp(zero_block, makeBlock(0, 1));

        int puncture_indx = 3;

        for (int c = 0; c < 3; c++)
        {
            // simulating opv sender
            prg.random_block(&seed, 1);
            auto t1 = clock_start();
            for (int i = 0; i < 100; i++)
            {
                opv_gen(seed, height, ot_msg_0, ot_msg_1, opv_sen);
            }
            std::cout << "\ntime from opv: " << time_from(t1) / 100 << " microseconds" << std::endl;

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
            get_choice_bits(puncture_indx, b, height - 2);
            for (int i = 0; i < height - 2; i++)
            {
                ot_rev_msg[i] = b[i] == 0 ? ot_msg_0[i] : ot_msg_1[i];
            }
            opv_reconstruction(b, ot_rev_msg, height, opv_rev);
            for (int i = 0; i < 2 * dimension; i++)
            {
                if (i != 2 * puncture_indx && i != (2 * puncture_indx + 1) && !cmpBlock(opv_sen + i, opv_rev + i, 1))
                {
                    std::cout << "ERROR - " << i << " - " << opv_sen[i] << " : " << opv_rev[i] << std::endl;
                }
            }
        }
        delete[] opv_sen;
        delete[] ot_msg_0;
        delete[] ot_msg_1;
        delete[] opv_rev;
        delete[] ot_rev_msg;
    }

    void test_OT(int length, int party, int port)
    {
        // NetIO *io = new NetIO(party==ALICE ? nullptr:"127.0.0.1", port);

        // IKNP<NetIO> * iknp = new IKNP<NetIO>(io);
        // std::cout <<"IKNP OT\t"<<double(length)/test_ot<IKNP<NetIO>>(iknp, io, party, length)*1e6<<" OTps"<<std::endl;

        block *msg_0, *msg_1, *rev_msg;
        ;
        msg_0 = new block[length];
        msg_1 = new block[length];
        rev_msg = new block[length];
        bool *b = new bool[length];

        this->prg.random_block(msg_0, length);
        this->prg.random_block(msg_1, length);
        this->prg.random_block(rev_msg, length);
        this->prg.random_bool(b, length);

        if (party == ALICE)
        {
            ot->send(msg_0, msg_1, length);
        }
        else
        {
            ot->recv(rev_msg, b, length);
        }
        io->flush();
        std::cout << "test_ot done.\n\n"
                  << std::endl;

        delete[] msg_0;
        delete[] msg_1;
        delete[] rev_msg;
        delete[] b;
        // delete iknp;
        // delete io;
    }

    void test_opv_with_OT()
    {
        // int length = (height-2);
        block *ot_sen_msg, *ot_rev_msg, *opv_sen, *opv_rev;
        opv_sen = new block[2 * dimension];
        opv_rev = new block[2 * dimension];
        ot_sen_msg = new block[2 * (height - 2)];
        ot_rev_msg = new block[height - 2];
        bool *b = new bool[height - 2];

        // FIXME: MUST INITILIZE BLOCKS, OTHERWISE YOU WILL GET ERRORS!
        this->prg.random_block(ot_sen_msg, 2 * (height - 2));
        this->prg.random_block(ot_rev_msg, (height - 2));
        // this->prg.random_bool(b, height-2);

        block seed = zero_block; // only for test
        // prg.random_block(&seed);

        if (party == ALICE)
        {
            this->opv_gen(seed, height, ot_sen_msg, ot_sen_msg + height - 2, opv_sen);
            ot->send(ot_sen_msg, ot_sen_msg + height - 2, height - 2);
        }
        else
        {
            int puncture_idx = 1;
            get_choice_bits(puncture_idx, b, height - 2);
            ot->recv(ot_rev_msg, b, height - 2);
            this->opv_reconstruction(b, ot_rev_msg, height, opv_rev);

            // do correctness check, only for test
            this->opv_gen(seed, height, ot_sen_msg, ot_sen_msg + height - 2, opv_sen);
            for (size_t i = 0; i < 2 * dimension; i++)
            {
                if (i != 2 * puncture_idx && i != (2 * puncture_idx + 1) && !cmpBlock(opv_sen + i, opv_rev + i, 1))
                {
                    std::cout << "ERROR - " << i << " - " << opv_sen[i] << " : " << opv_rev[i] << std::endl;
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

    void test_opm()
    {

        auto t1 = clock_start();
        this->opm_gen(height, dimension, ggm_master_key, opm_sen_memory, ot_sen_msg, ot_sen_msg + (height - 2) * dimension, opm_a, opm_b, column_sum);
        std::cout << "time from opm_gen(): " << time_from(t1) << " microseconds" << std::endl;

        // simulating OTs
        block *ot_sen_msg_0 = this->ot_sen_msg, *ot_sen_msg_1 = this->ot_sen_msg + (height - 2) * dimension;
        bool b[height - 2];
        for (size_t i = 0; i < dimension; i++)
        {
            get_choice_bits(permutation[i], b, height - 2);
            for (int j = 0; j < height - 2; j++)
            {
                ot_rev_msg[i * (height - 2) + j] = b[j] == false ? ot_sen_msg_0[i * (height - 2) + j] : ot_sen_msg_1[i * (height - 2) + j]; // NOTE: b[j] not b[i] !!!
            }
        }

        auto t2 = clock_start();
        this->opm_reconstruction(height, dimension, permutation, permutation_inverse, ot_rev_msg, column_sum, opm_rev_memory, opm_delta);
        std::cout << "time from opm_reconstruction(): " << time_from(t2) << " microseconds" << std::endl;

        // simulating consistency check
        block *chi = new block[dimension * dimension];
        Hash hash;
        block digest[2];

        auto t3 = clock_start();
        // consistency check method 1
        hash.hash_once(digest, this->column_sum, dimension * sizeof(block));
        uni_hash_coeff_gen(chi, digest[0], dimension * dimension);
        for (size_t i = 0; i < dimension; i++)
        {
            std::cout << " column sum i = " << i << " : " << column_sum[i];
        }
        std::cout << std::endl;
        std::cout << "digest[0]: " << digest[0] << std::endl;
        std::cout << "chi[0]: " << chi[0] << std::endl;

        block r_sen = zero_block, r_rev = zero_block;
        block r1, r2;
        for (size_t i = 0; i < dimension; i++)
        {
            for (size_t j = 0; j < dimension; j++)
            {
                gfmul(this->opm_sen_memory[2 * i * dimension + 2 * j], chi[i * dimension + j], &r1);
                gfmul(this->opm_rev_memory[2 * i * dimension + 2 * j], chi[i * dimension + j], &r2);
                r_sen = r_sen ^ r1;
                r_rev = r_rev ^ r2;
            }
        }
        if (!cmpBlock(&r_sen, &r_rev, 1))
        {
            log::error("OPM check fails!\n");
        }
        std::cout << "r_sen: " << r_sen << ", r_rev: " << r_rev << std::endl;
        std::cout << "time from OPM check method 1 : " << time_from(t3) << " microseconds" << std::endl;

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
        for (size_t i = 0; i < dimension; i++)
        {
            // delta = p(a)+b
            opt_check[i] = opm_delta[i] ^ opm_a[permutation[i]] ^ opm_b[i];
            if (!cmpBlock(opt_check + i, &zero_block, 1))
            {
                log::error("opt check fails.\n");
            }
        }
        delete[] opt_check;

        delete[] chi;
    }

    void test_opm_with_OT()
    {
        block *chi = new block[dimension * dimension];
        Hash hash;
        block digest[2];
        block sen_tag = zero_block, rev_tag = zero_block;

        if (party == ALICE)
        {
            // auto t1 = clock_start();
            this->opm_gen(height, dimension, ggm_master_key, opm_sen_memory, ot_sen_msg, ot_sen_msg + (height - 2) * dimension, opm_a, opm_b, column_sum);
            // std::cout << "time from opm_gen(): " << time_from(t1) << " microseconds" <<std::endl;

            hash.hash_once(digest, column_sum, dimension * sizeof(block)); //`column_sum` was updated by function `opm_gen()`
            uni_hash_coeff_gen(chi, digest[0], dimension * dimension);
            block r1;
            for (size_t i = 0; i < dimension; i++)
            {
                for (size_t j = 0; j < dimension; j++)
                {
                    gfmul(this->opm_sen_memory[2 * i * dimension + 2 * j], chi[i * dimension + j], &r1);
                    sen_tag = sen_tag ^ r1;
                }
            }

            ot->send(ot_sen_msg, ot_sen_msg + (height - 2) * dimension, (height - 2) * dimension);
            // std::cout<< "ot sen msg: " << ot_sen_msg[0] << " " << ot_sen_msg[(height-2)*dimension] <<std::endl;
            io->send_data(column_sum, dimension * sizeof(block));
            io->flush();
            // std::cout<<"sen column_sum: "<< *column_sum <<std::endl;
            hash.hash_once(digest, &sen_tag, sizeof(block));
            io->send_data(digest, 2 * sizeof(block));
            io->flush();
            // std::cout<<"sen digest: "<< digest[0] << " " << digest[1] <<std::endl;
            io->recv_data(&rev_tag, sizeof(block)); // io->flush();
            // do check
            if (!cmpBlock(&rev_tag, &sen_tag, 1))
            {
                log::error("hash from receiver does not match sender's tag\n");
            }
        }
        else
        {
            for (size_t i = 0; i < dimension; i++)
            {
                get_choice_bits(permutation[i], ot_choise_bits + i * (height - 2), height - 2);
            }

            ot->recv(ot_rev_msg, ot_choise_bits, (height - 2) * dimension);
            // std::cout<< "ot rev msg: " << *ot_rev_msg <<std::endl;
            io->recv_data(column_sum, dimension * sizeof(block));

            // io->flush();

            // auto t2 = clock_start();
            this->opm_reconstruction(height, dimension, permutation, permutation_inverse, ot_rev_msg, column_sum, opm_rev_memory, opm_delta);
            // std::cout << "time from opm_reconstruction(): " << time_from(t2) << " microseconds" <<std::endl;
            // std::cout<<"rev column_sum: "<< *column_sum <<std::endl;
            //  consistency check
            Hash hash;
            block digest[2];

            // auto t3 = clock_start();
            hash.hash_once(digest, column_sum, dimension * sizeof(block));
            uni_hash_coeff_gen(chi, digest[0], dimension * dimension);
            block r1;
            for (size_t i = 0; i < dimension; i++)
            {
                for (size_t j = 0; j < dimension; j++)
                {
                    gfmul(this->opm_rev_memory[2 * i * dimension + 2 * j], chi[i * dimension + j], &r1);
                    rev_tag = rev_tag ^ r1;
                }
            }

            block hash_sen[2];
            io->recv_data(hash_sen, 2 * sizeof(block));
            hash.hash_once(digest, &rev_tag, sizeof(block));
            // std::cout<<"rev digest: "<< digest[0] << " " << digest[1] <<std::endl;
            if (!cmpBlock(hash_sen, digest, 2))
            {
                log::error("hash from sender does not match receiver's tag\n");
                // return ;//FIXME:
            }

            // send back the tag for check
            io->send_data(&rev_tag, sizeof(block));
            io->flush();

            // std::cout << "time from OPM check method 1 : " << time_from(t3) << " microseconds" <<std::endl;
        }
        // io->flush();

        delete[] chi;
        std::cout << "test_opm_with_OT DONE." << std::endl;
    }

    int get_gbn_batch_size(const size_t dimension, const size_t width)
    {
        int n_bits_dimension, n_bits_width, d;
        n_bits_dimension = log2(dimension);
        n_bits_width = log2(width);
        // std::cout << (1.0 * n_bits_width / n_bits_dimension) << std::endl;
        d = int(ceil(1.0 * n_bits_width / n_bits_dimension));

        return (2 * d - 1) * (width / dimension);
    }

    void test_batch_semi_opm_with_OT(const size_t width, const size_t demension, int batch_size)
    {
        // update batch_size
        // batch_size = get_gbn_batch_size(dimension, width);

        std::cout << "bits_dimension = " << int(log2(dimension)) << ", bits_width = " << int(log2(width)) << ", batch OPTs size: " << batch_size << std::endl;

        // sender memory
        block *batch_ggm_master_key, *batch_ot_sen_msg, *batch_opm_a, *batch_opm_b;
        block *opm_sen_memory;

        // receiver memory
        block *batch_ot_rev_msg, *batch_opm_delta, *opm_rev_memory;
        std::vector<std::vector<size_t>> batch_permutations, batch_inverse_permutations;
        bool *batch_ot_choise_bits;

        auto t1 = clock_start();
        if (party == ALICE)
        { // for sender-side batch memory

            batch_ggm_master_key = new block[batch_size * dimension];
            prg.random_block(batch_ggm_master_key, batch_size * dimension);
            // batch_opm_sen_memory = new block[batch_size * 2 * dimension * dimension];
            batch_ot_sen_msg = new block[batch_size * 2 * (height - 1) * dimension];
            batch_opm_a = new block[batch_size * dimension];
            batch_opm_b = new block[batch_size * dimension];

            // int total_size = ;
            std::cout << batch_size * dimension * dimension << std::endl;

            // for sender-side instance memory
            block *ggm_master_key;
            block *ot_sen_msg_0, *ot_sen_msg_1;
            block *opm_a;
            block *opm_b;

            // run local opm generation
            std::cout << "begin opm_gen()" << std::endl;
            auto t_opm = clock_start();

            opm_sen_memory = new block[dimension * dimension];
            Hash hash;

            for (size_t i = 0; i < batch_size; i++)
            {
                ggm_master_key = &batch_ggm_master_key[i * dimension];
                // opm_sen_memory = &batch_opm_sen_memory[i * 2 * dimension * dimension];
                ot_sen_msg_0 = &batch_ot_sen_msg[i * (height - 1) * dimension];
                ot_sen_msg_1 = ot_sen_msg_0 + batch_size * (height - 1) * dimension;
                opm_a = &batch_opm_a[i * dimension];
                opm_b = &batch_opm_b[i * dimension];
                // opm_gen
                for (size_t row = 0; row < dimension; row++)
                {
                    this->ggm_tree_gen(ggm_master_key[row], height, ot_sen_msg_0, ot_sen_msg_1, &opm_sen_memory[row * dimension]);
                }

                // set opm_a, opm_b
                for (size_t row = 0; row < dimension; row++)
                {
                    for (int col = 0; col < dimension; col++)
                    {
                        opm_b[row] = opm_b[row] ^ opm_sen_memory[row * dimension + col];
                        opm_a[col] = opm_a[col] ^ opm_sen_memory[row * dimension + col];
                    }
                }
            } // end for

            // std::cout<< batch_size * 2 * dimension * dimension <<std::endl;
            std::cout << "semi opm_gen() takes " << time_from(t_opm) << std::endl;

            // do ot and io
            std::cout << "begin ot and io" << std::endl;
            auto t_ot = clock_start();
            // send ot msgs
            ot->send(batch_ot_sen_msg, batch_ot_sen_msg + batch_size * (height - 1) * dimension, batch_size * (height - 1) * dimension);
            std::cout << "semi ot and io takes " << time_from(t_ot) << std::endl;
        } // end sender
        else
        {
            batch_ot_rev_msg = new block[batch_size * (height - 1) * dimension];
            batch_opm_delta = new block[batch_size * dimension];
            batch_ot_choise_bits = new bool[batch_size * (height - 1) * dimension];

            // prepare receiver permutations and ot choice bits
            std::cout << "begin permutation prepare" << std::endl;
            auto t_p = clock_start();
            for (size_t i = 0; i < batch_size; i++)
            {
                bool *ot_choise_bits = &batch_ot_choise_bits[i * (height - 1) * dimension];

                std::vector<size_t> permutation, permutation_inverse;
                get_random_permutation(dimension, permutation, permutation_inverse);
                batch_permutations.push_back(permutation);
                batch_inverse_permutations.push_back(permutation_inverse);

                for (size_t row = 0; row < dimension; row++)
                {
                    get_choice_bits(permutation[row], ot_choise_bits + row * (height - 1), height - 1);
                }
            }
            std::cout << "semi pre permutaiton takes " << time_from(t_p) << std::endl;

            std::cout << "begin ot io" << std::endl;
            auto t_ot = clock_start();
            // receive ot msgs and column sum
            ot->recv(batch_ot_rev_msg, batch_ot_choise_bits, batch_size * (height - 1) * dimension);
            // std::cout<< "malicious bob here 1" <<std::endl;
            std::cout << "ot and io taks " << time_from(t_ot) << std::endl;

            // recover batch OPMs
            std::cout << "begin opm rec()" << std::endl;
            auto t_rec = clock_start();
            opm_rev_memory = new block[dimension * dimension];
            for (size_t i = 0; i < batch_size; i++)
            {
                block *ot_rev_msg = &batch_ot_rev_msg[i * (height - 1) * dimension];
                bool *ot_choise_bits = &batch_ot_choise_bits[i * (height - 1) * dimension];
                for (size_t row = 0; row < dimension; row++)
                {
                    this->ggm_tree_reconstruction(&ot_choise_bits[row * (height - 1)], ot_rev_msg, height, &opm_rev_memory[row * dimension]);
                }
                // set opm_delta
                block *opm_delta = &batch_opm_delta[i * dimension];
                for (size_t row = 0; row < dimension; row++)
                {
                    for (size_t col = 0; col < dimension; col++)
                    {
                        if (col != permutation[row])
                        {
                            opm_delta[row] ^= opm_rev_memory[row * dimension + col];
                            opm_delta[permutation_inverse[col]] ^= opm_rev_memory[row * dimension + col];
                        }
                    }
                }
            }
            std::cout << "opm_reconstruction taks " << time_from(t_rec) << std::endl;
        } // end receiver

        std::cout << "bits_dimension = " << int(log2(dimension)) << ", bits_width = " << int(log2(width)) << ". Time for OPT generation: " << time_from(t1) << " millseconds" << std::endl;
        std::cout << "total communication size: " << io->counter << " bytes" << std::endl;

        //---------------------------------------write to file-------------------------------------------
        // std::cout << "opm file store begins." << std::endl;

        // std::string filename = party == ALICE ? OPT_SENDER_FILE : OPT_RECEIVER_FILE;
        // filename = filename + "_OFF_SEMI_TEST_" + std::to_string(height - 1) + "_" + std::to_string(int(log2(width)));
        // std::cout << "filename: " << filename << std::endl;

        // auto t_file = clock_start();
        // std::ofstream outfile(filename);
        // if (outfile.is_open())
        //     outfile.close();
        // else
        //     error("create a directory to store opm data");
        // FileIO fio(filename.c_str(), false);
        // fio.send_data(&party, sizeof(int));
        // fio.send_data(&batch_size, sizeof(int));

        // for (size_t i = 0; i < batch_size; i++)
        // {
        //     if (party == ALICE)
        //     { // ALICE is the OPM sender
        //         fio.send_data(&batch_opm_a[i * dimension], sizeof(block) * dimension);
        //         fio.send_data(&batch_opm_b[i * dimension], sizeof(block) * dimension);
        //     }
        //     else
        //     {
        //         fio.send_data(batch_permutations[i].data(), sizeof(size_t) * dimension);
        //         fio.send_data(&batch_opm_delta[i * dimension], sizeof(block) * dimension);
        //     }
        // }
        // std::cout << "time for file store: " << time_from(t_file) << std::endl;

        if (party == ALICE)
        {
            delete[] batch_ggm_master_key;
            delete[] batch_opm_a;
            delete[] batch_opm_b;
            delete[] batch_ot_sen_msg;
        }
        else{
            delete[] opm_rev_memory;
            delete[] batch_ot_rev_msg;
            delete[] batch_opm_delta;
            delete[] batch_ot_choise_bits;
        }
    }

    void test_batch_opm_with_OT(const size_t width, const size_t demension, int batch_size)
    {
        // update batch_size
        // batch_size = get_gbn_batch_size(dimension, width);

        std::cout << "bits_dimension = " << int(log2(dimension)) << ", bits_width = " << int(log2(width)) << ", batch OPTs size: " << batch_size << std::endl;

        // sender memory
        block *batch_ggm_master_key, *batch_ot_sen_msg, *batch_opm_a, *batch_opm_b;

        // receiver memory
        block *batch_ot_rev_msg, *batch_opm_delta;
        std::vector<std::vector<size_t>> batch_permutations, batch_inverse_permutations;
        bool *batch_ot_choise_bits;

        // common for sender and receiver
        block *batch_column_sum = new block[batch_size * dimension];

        auto t1 = clock_start();
        if (party == ALICE)
        { // for sender-side batch memory

            batch_ggm_master_key = new block[batch_size * dimension];
            // batch_opm_sen_memory = new block[batch_size * 2 * dimension * dimension];
            batch_ot_sen_msg = new block[batch_size * 2 * (height - 2) * dimension];
            batch_opm_a = new block[batch_size * dimension];
            batch_opm_b = new block[batch_size * dimension];

            // int total_size = ;
            std::cout << batch_size * 2 * dimension * dimension << std::endl;

            // for sender-side instance memory
            block *ggm_master_key;
            block *opm_sen_memory;
            block *ot_sen_msg;
            block *opm_a;
            block *opm_b;
            block *column_sum;

            // run local opm generation
            std::cout << "begin opm_gen()" << std::endl;
            auto t_opm = clock_start();

            opm_sen_memory = new block[2 * dimension * dimension];
            Hash hash;

#ifdef MALICIOUS
            block *chi = new block[dimension * dimension];
            block digest[2];
            block sen_tag = zero_block;
#endif

            for (size_t i = 0; i < batch_size; i++)
            {
                ggm_master_key = &batch_ggm_master_key[i * dimension];
                // opm_sen_memory = &batch_opm_sen_memory[i * 2 * dimension * dimension];
                ot_sen_msg = &batch_ot_sen_msg[i * (height - 2) * dimension];
                opm_a = &batch_opm_a[i * dimension];
                opm_b = &batch_opm_b[i * dimension];
                column_sum = &batch_column_sum[i * dimension];
                // opm_gen
                this->opm_gen(height, dimension, ggm_master_key, opm_sen_memory, ot_sen_msg, ot_sen_msg + batch_size * (height - 2) * dimension, opm_a, opm_b, column_sum);

// updtae opm tag
#ifdef MALICIOUS
                hash.hash_once(digest, column_sum, dimension * sizeof(block)); //`column_sum` was updated by function `opm_gen()`
                uni_hash_coeff_gen(chi, digest[0], dimension * dimension);
                // compute sen_tag
                block r1;
                for (size_t j = 0; j < dimension; j++)
                {
                    for (size_t k = 0; k < dimension; k++)
                    {
                        gfmul(opm_sen_memory[2 * j * dimension + 2 * k], chi[j * dimension + k], &r1);
                        sen_tag = sen_tag ^ r1;
                    }
                }

#endif
            } // end for

#ifdef MALICIOUS
            delete[] chi;
#endif
            delete[] opm_sen_memory;
            delete[] batch_ggm_master_key;

            // std::cout<< batch_size * 2 * dimension * dimension <<std::endl;
            std::cout << "opm_gen() takes " << time_from(t_opm) << std::endl;

            // do ot and io
            std::cout << "begin ot and io" << std::endl;
            auto t_ot = clock_start();
            // send ot msgs
            ot->send(batch_ot_sen_msg, batch_ot_sen_msg + batch_size * (height - 2) * dimension, batch_size * (height - 2) * dimension);
            std::cout << batch_size * (height - 2) * dimension << " : " << "ot takes " << time_from(t_ot) << std::endl;

            // send column-wise sum
            auto t_io = clock_start();
            io->send_data(batch_column_sum, batch_size * dimension * sizeof(block));
            io->flush();
            std::cout << batch_size * dimension << " : " << "io takes " << time_from(t_io) << std::endl;

#ifdef MALICIOUS
            // compare equaity tag
            block sen_tag_hash[2];
            hash.hash_once(sen_tag_hash, &sen_tag, sizeof(block));
            io->send_data(sen_tag_hash, 2 * sizeof(block));
            io->flush();
            // std::cout<< "malicious alice here 1" <<std::endl;

            block rev_tag = zero_block;
            io->recv_data(&rev_tag, sizeof(block)); // io->flush();

            // std::cout << "malicious alice here2" << std::endl;
            //  do check
            if (!cmpBlock(&rev_tag, &sen_tag, 1))
            {
                log::error("hash from receiver does not match sender's tag\n");
            }

            std::cout << "opm gen after check done." << std::endl;
#endif
        }
        else
        { // for receiver-side batch memory
            batch_ot_rev_msg = new block[batch_size * (height - 2) * dimension];
            // batch_opm_rev_memory = new block[batch_size * 2 * dimension * dimension];
            batch_opm_delta = new block[batch_size * dimension];
            batch_ot_choise_bits = new bool[batch_size * (height - 2) * dimension];

            // for receiver-side instance memory
            block *ot_rev_msg;
            block *column_sum;
            block *opm_rev_memory;
            block *opm_delta;

            // prepare receiver permutations and ot choice bits
            std::cout << "begin permutation prepare" << std::endl;
            auto t_p = clock_start();
            for (size_t i = 0; i < batch_size; i++)
            {
                bool *ot_choise_bits = &batch_ot_choise_bits[i * (height - 2) * dimension];

                std::vector<size_t> permutation, permutation_inverse;
                get_random_permutation(dimension, permutation, permutation_inverse);
                batch_permutations.push_back(permutation);
                batch_inverse_permutations.push_back(permutation_inverse);

                for (size_t j = 0; j < dimension; j++)
                {
                    get_choice_bits(permutation[j], ot_choise_bits + j * (height - 2), height - 2);
                }
                // std::cout<< "malicious bob here 0" <<std::endl;
            }
            std::cout << "pre permutaiton takes " << time_from(t_p) << std::endl;

            std::cout << "begin ot io" << std::endl;
            auto t_ot = clock_start();
            // receive ot msgs and column sum
            ot->recv(batch_ot_rev_msg, batch_ot_choise_bits, batch_size * (height - 2) * dimension);
            // std::cout<< "malicious bob here 1" <<std::endl;
            std::cout << batch_size * (height - 2) * dimension << " : " << "ot takes " << time_from(t_ot) << std::endl;

            auto t_io = clock_start();
            io->recv_data(batch_column_sum, batch_size * dimension * sizeof(block));
            io->flush();
            // std::cout<< "malicious bob here 2" <<std::endl;
            std::cout << batch_size * dimension << " : " << "ot and io taks " << time_from(t_io) << std::endl;

            // recover batch OPMs
            std::cout << "begin opm rec()" << std::endl;
            auto t_rec = clock_start();

            opm_rev_memory = new block[2 * dimension * dimension];
#ifdef MALICIOUS
            block *chi = new block[dimension * dimension];
            Hash hash;
            block rev_tag = zero_block;
            block digest[2];
#endif

            for (size_t i = 0; i < batch_size; i++)
            {
                ot_rev_msg = &batch_ot_rev_msg[i * (height - 2) * dimension];
                column_sum = &batch_column_sum[i * dimension];
                opm_delta = &batch_opm_delta[i * dimension];

                this->opm_reconstruction(height, dimension, batch_permutations[i], batch_inverse_permutations[i], ot_rev_msg, column_sum, opm_rev_memory, opm_delta);
                // std::cout<< "malicious bob here 3" <<std::endl;

#ifdef MALICIOUS

                hash.hash_once(digest, column_sum, dimension * sizeof(block));
                uni_hash_coeff_gen(chi, digest[0], dimension * dimension);
                // std::cout<< "malicious bob here 4.0" <<std::endl;
                block r1;
                for (size_t j = 0; j < dimension; j++)
                {
                    for (size_t k = 0; k < dimension; k++)
                    {
                        gfmul(opm_rev_memory[2 * j * dimension + 2 * k], chi[j * dimension + k], &r1);
                        rev_tag = rev_tag ^ r1;
                    }
                }
#endif
            } // end for
            std::cout << "opm_reconstruction taks " << time_from(t_rec) << std::endl;
            delete[] opm_rev_memory;

#ifdef MALICIOUS
            block sen_tag_hash[2];
            io->recv_data(sen_tag_hash, 2 * sizeof(block));
            // std::cout<< "malicious bob here 5" <<std::endl;

            block rev_tag_hash[2];
            hash.hash_once(rev_tag_hash, &rev_tag, sizeof(block));
            // std::cout<<"rev digest: "<< digest[0] << " " << digest[1] <<std::endl;
            if (!cmpBlock(sen_tag_hash, rev_tag_hash, 2))
            {
                log::error("hash from sender does not match receiver's tag\n");
                // return ;//FIXME:
            }

            // send back rev_tag for check
            io->send_data(&rev_tag, sizeof(block));
            io->flush();
#endif

#ifdef MALICIOUS
            delete[] chi;
#endif
        } // end BOB

        std::cout << "bits_dimension = " << int(log2(dimension)) << ", bits_width = " << int(log2(width)) << ". Time for OPT generation: " << time_from(t1) << " millseconds" << std::endl;
        std::cout << "total communication size: " << io->counter << " bytes" << std::endl;
        //---------------------------------------write to file-------------------------------------------
        std::cout << "opm file store begins." << std::endl;

        std::string filename = party == ALICE ? OPT_SENDER_FILE : OPT_RECEIVER_FILE;
        filename = filename + "_OFF_MALI_TEST_" + std::to_string(height - 2) + "_" + std::to_string(int(log2(width)));
        std::cout << "filename: " << filename << std::endl;

        auto t_file = clock_start();
        std::ofstream outfile(filename);
        if (outfile.is_open())
            outfile.close();
        else
            error("create a directory to store opm data");
        FileIO fio(filename.c_str(), false);
        fio.send_data(&party, sizeof(int));
        fio.send_data(&batch_size, sizeof(int));

        for (size_t i = 0; i < batch_size; i++)
        {
            if (party == ALICE)
            { // ALICE is the OPM sender
                fio.send_data(&batch_opm_a[i * dimension], sizeof(block) * dimension);
                fio.send_data(&batch_opm_b[i * dimension], sizeof(block) * dimension);
            }
            else
            {
                fio.send_data(batch_permutations[i].data(), sizeof(size_t) * dimension);
                fio.send_data(&batch_opm_delta[i * dimension], sizeof(block) * dimension);
            }
        }
        std::cout << "time for file store: " << time_from(t_file) << std::endl;

        //---------------------------------------------------------------------------------------
        // clean-up
        if (party == ALICE)
        { // for sender-side memory
            // delete[] batch_ggm_master_key;
            delete[] batch_opm_a;
            delete[] batch_opm_b;
        }
        else
        { // for receiver-side memory
            delete[] batch_ot_rev_msg;
            delete[] batch_opm_delta;
        }
        delete[] batch_column_sum;
    }

    void test_batch_opm_with_OT_threads(const size_t width, const int height, const size_t dimension, const bool malicious = true, int batch_size = MAX_BATCH_OPM_SIZE)
    {
        // update batch_size
        batch_size = get_gbn_batch_size(dimension, width);
        std::cout << "batch OPTs size: " << batch_size << std::endl;

        // sender memory
        block *batch_ggm_master_key, *batch_opm_sen_memory, *batch_ot_sen_msg, *batch_opm_a, *batch_opm_b;

        // receiver memory
        block *batch_ot_rev_msg, *batch_opm_rev_memory, *batch_opm_delta;
        std::vector<std::vector<size_t>> batch_permutations, batch_inverse_permutations;
        bool *batch_ot_choise_bits;

        // common for sender and receiver
        block *batch_column_sum = new block[batch_size * dimension];

        auto t1 = clock_start();
        if (party == ALICE)
        { // for sender-side batch memory

            batch_ggm_master_key = new block[batch_size * dimension];
            batch_opm_sen_memory = new block[batch_size * 2 * dimension * dimension];
            batch_ot_sen_msg = new block[batch_size * 2 * (height - 2) * dimension];
            batch_opm_a = new block[batch_size * dimension];
            batch_opm_b = new block[batch_size * dimension];

            // for sender-side instance memory

            // run local opm generation
            std::vector<std::future<void>> fut;
            int sub_batch_size = batch_size / n_threads;
            for (size_t t = 0; t < n_threads; t++)
            {
                size_t start, end;
                start = t * sub_batch_size;
                end = (t + 1) * sub_batch_size;
                if (t = n_threads - 1)
                    end = batch_size;

                fut.push_back(this->pool->enqueue(
                    [this, height, dimension, batch_size, t, start, end, batch_ggm_master_key, batch_opm_sen_memory, batch_ot_sen_msg, batch_opm_a, batch_opm_b, batch_column_sum]()
                    {
                        for (size_t i = start; i < end; i++)
                        {
                            block *ggm_master_key = &batch_ggm_master_key[i * dimension];
                            block *opm_sen_memory = &batch_opm_sen_memory[i * 2 * dimension * dimension];
                            block *ot_sen_msg = &batch_ot_sen_msg[i * (height - 2) * dimension];
                            block *opm_a = &batch_opm_a[i * dimension];
                            block *opm_b = &batch_opm_b[i * dimension];
                            block *column_sum = &batch_column_sum[i * dimension];
                            // opm_gen
                            this->opm_gen(height, dimension, ggm_master_key, opm_sen_memory, ot_sen_msg, ot_sen_msg + batch_size * (height - 2) * dimension, opm_a, opm_b, column_sum);
                        }
                    }));
            }
            for (auto &f : fut)
                f.get();

            // send ot msgs
            ot->send(batch_ot_sen_msg, batch_ot_sen_msg + batch_size * (height - 2) * dimension, batch_size * (height - 2) * dimension);

            // send column-wise sum
            io->send_data(batch_column_sum, batch_size * dimension * sizeof(block));
            io->flush();

            // do consistency check
            if (malicious)
            {
                std::vector<std::future<void>> fut;
                int sub_batch_size = batch_size / n_threads;
                Hash hash;
                block sen_tags[n_threads];
                memset(sen_tags, 0, n_threads * sizeof(block));

                for (size_t t = 0; t < n_threads; t++)
                {
                    size_t start, end;
                    start = t * sub_batch_size;
                    end = (t + 1) * sub_batch_size;
                    if (t = n_threads - 1)
                        end = batch_size;

                    fut.push_back(this->pool->enqueue(
                        [this, height, dimension, batch_size, t, start, end, batch_opm_sen_memory, batch_column_sum, &hash, &sen_tags]()
                        {
                            block *chi = new block[dimension * dimension];
                            block digest[2];
                            // block *sub_sen_tag = sen_tags + t;

                            for (size_t i = start; i < end; i++)
                            {
                                block *opm_sen_memory = &batch_opm_sen_memory[i * 2 * dimension * dimension];
                                block *column_sum = &batch_column_sum[i * dimension];

                                // opm_check
                                hash.hash_once(digest, column_sum, dimension * sizeof(block)); //`column_sum` was updated by function `opm_gen()`
                                uni_hash_coeff_gen(chi, digest[0], dimension * dimension);

                                // compute sen_tag
                                block r1;
                                for (size_t j = 0; j < dimension; j++)
                                {
                                    for (size_t k = 0; k < dimension; k++)
                                    {
                                        gfmul(opm_sen_memory[2 * j * dimension + 2 * k], chi[j * dimension + k], &r1);
                                        sen_tags[t] = sen_tags[t] ^ r1;
                                    }
                                }
                            }
                        }));
                }
                for (auto &f : fut)
                    f.get();
                for (int t = 1; t < n_threads; t++)
                {
                    sen_tags[0] = sen_tags[0] ^ sen_tags[t];
                }

                // compare equaity tag
                block sen_tag_hash[2];
                hash.hash_once(sen_tag_hash, &sen_tags[0], sizeof(block));
                io->send_data(sen_tag_hash, 2 * sizeof(block));
                io->flush();
                // std::cout<< "malicious alice here 1" <<std::endl;

                block rev_tag;
                io->recv_data(&rev_tag, sizeof(block)); // io->flush();

                std::cout << "malicious alice here2" << std::endl;
                // do check
                if (!cmpBlock(&rev_tag, &sen_tags[0], 1))
                {
                    log::error("hash from receiver does not match sender's tag\n");
                }
                // std::cout<< "malicious alice here 3" <<std::endl;
            }
        }
        else
        { // for receiver-side batch memory
            batch_ot_rev_msg = new block[batch_size * (height - 2) * dimension];
            batch_opm_rev_memory = new block[batch_size * 2 * dimension * dimension];
            batch_opm_delta = new block[batch_size * dimension];
            batch_ot_choise_bits = new bool[batch_size * (height - 2) * dimension];

            // for receiver-side instance memory
            block *ot_rev_msg;
            block *column_sum;
            block *opm_rev_memory;
            block *opm_delta;

            // prepare receiver permutations and ot choice bits
            for (size_t i = 0; i < batch_size; i++)
            {
                bool *ot_choise_bits = &batch_ot_choise_bits[i * (height - 2) * dimension];

                std::vector<size_t> permutation, permutation_inverse;
                get_random_permutation(dimension, permutation, permutation_inverse);
                batch_permutations.push_back(permutation);
                batch_inverse_permutations.push_back(permutation_inverse);

                for (size_t j = 0; j < dimension; j++)
                {
                    get_choice_bits(permutation[j], ot_choise_bits + j * (height - 2), height - 2);
                }
                // std::cout<< "malicious bob here 0" <<std::endl;
            }

            // receive ot msgs and column sum
            ot->recv(batch_ot_rev_msg, batch_ot_choise_bits, batch_size * (height - 2) * dimension);
            // std::cout<< "malicious bob here 1" <<std::endl;

            io->recv_data(batch_column_sum, batch_size * dimension * sizeof(block));
            io->flush();
            // std::cout<< "malicious bob here 2" <<std::endl;

            // recover batch OPMs
            std::vector<std::future<void>> fut;
            size_t sub_batch_size = batch_size / n_threads;
            for (size_t t = 0; t < n_threads; t++)
            {
                size_t start, end;
                start = t * sub_batch_size;
                end = (t + 1) * sub_batch_size;
                if (t = n_threads - 1)
                    end = batch_size;

                fut.push_back(this->pool->enqueue(
                    [this, height, dimension, t, start, end, batch_ot_rev_msg, batch_column_sum, batch_opm_rev_memory, batch_opm_delta, &batch_permutations, &batch_inverse_permutations]()
                    {
                        for (size_t i = start; i < end; i++)
                        {
                            block *ot_rev_msg = &batch_ot_rev_msg[i * (height - 2) * dimension];
                            block *column_sum = &batch_column_sum[i * dimension];
                            block *opm_rev_memory = &batch_opm_rev_memory[i * 2 * dimension * dimension];
                            block *opm_delta = &batch_opm_delta[i * dimension];

                            this->opm_reconstruction(height, dimension, batch_permutations[i], batch_inverse_permutations[i], ot_rev_msg, column_sum, opm_rev_memory, opm_delta);
                        }
                    }));
            }
            for (auto &f : fut)
                f.get();

            // do consistency check
            if (malicious)
            {
                std::vector<std::future<void>> fut;
                int sub_batch_size = batch_size / n_threads;
                Hash hash;
                block rev_tags[n_threads];
                memset(rev_tags, 0, n_threads * sizeof(block));

                for (size_t t = 0; t < n_threads; t++)
                {
                    size_t start, end;
                    start = t * sub_batch_size;
                    end = (t + 1) * sub_batch_size;
                    if (t = n_threads - 1)
                        end = batch_size;

                    fut.push_back(this->pool->enqueue(
                        [this, height, dimension, batch_size, t, start, end, batch_opm_rev_memory, batch_column_sum, &hash, &rev_tags]()
                        {
                            block *chi = new block[dimension * dimension];
                            block digest[2];
                            // block *sub_sen_tag = sen_tags + t;

                            for (size_t i = start; i < end; i++)
                            {
                                block *opm_rev_memory = &batch_opm_rev_memory[i * 2 * dimension * dimension];
                                block *column_sum = &batch_column_sum[i * dimension];

                                // opm_check
                                hash.hash_once(digest, column_sum, dimension * sizeof(block)); //`column_sum` was updated by function `opm_gen()`
                                uni_hash_coeff_gen(chi, digest[0], dimension * dimension);

                                // compute sen_tag
                                block r1;
                                for (size_t j = 0; j < dimension; j++)
                                {
                                    for (size_t k = 0; k < dimension; k++)
                                    {
                                        gfmul(opm_rev_memory[2 * j * dimension + 2 * k], chi[j * dimension + k], &r1);
                                        rev_tags[t] = rev_tags[t] ^ r1;
                                    }
                                }
                            }
                        }));
                }
                for (auto &f : fut)
                    f.get();
                for (size_t t = 1; t < n_threads; t++)
                {
                    rev_tags[0] = rev_tags[0] ^ rev_tags[t];
                }

                block sen_tag_hash[2], rev_tag_hash[2];
                io->recv_data(sen_tag_hash, 2 * sizeof(block));
                // std::cout<< "malicious bob here 5" <<std::endl;
                hash.hash_once(rev_tag_hash, &rev_tags[0], sizeof(block));
                // std::cout<<"rev digest: "<< digest[0] << " " << digest[1] <<std::endl;
                if (!cmpBlock(sen_tag_hash, rev_tag_hash, 2))
                {
                    log::error("hash from sender does not match receiver's tag\n");
                    // return ;//FIXME:
                }

                // send back rev_tag for check
                io->send_data(&rev_tags[0], sizeof(block));
                io->flush();
                // std::cout<< "malicious bob here 6" <<std::endl;
            }
        }
        std::cout << "bits_dimension = " << int(log2(dimension)) << ", bits_width = " << int(log2(width)) << ". Time for OPT generation: " << time_from(t1) << std::endl;

        //---------------------------------------write to file-------------------------------------------
        std::string filename = party == ALICE ? OPT_SENDER_FILE : OPT_RECEIVER_FILE;
        filename = filename + "_OFF_TEST_" + std::to_string(height - 2) + "_" + std::to_string(int(log2(width)));
        std::ofstream outfile(filename);
        if (outfile.is_open())
            outfile.close();
        else
            error("create a directory to store opm data");
        FileIO fio(filename.c_str(), false);
        fio.send_data(&party, sizeof(int));
        fio.send_data(&batch_size, sizeof(int));

        for (size_t i = 0; i < batch_size; i++)
        {
            if (party == ALICE)
            { // ALICE is the OPM sender
                fio.send_data(&batch_opm_a[i * dimension], sizeof(block) * dimension);
                fio.send_data(&batch_opm_b[i * dimension], sizeof(block) * dimension);
            }
            else
            {
                fio.send_data(batch_permutations[i].data(), sizeof(size_t) * dimension);
                fio.send_data(&batch_opm_delta[i * dimension], sizeof(block) * dimension);
            }
        }
        //---------------------------------------------------------------------------------------
        // clean-up
        if (party == ALICE)
        { // for sender-side memory
            delete[] batch_ggm_master_key;
            delete[] batch_opm_sen_memory;
            delete[] batch_ot_sen_msg;
            delete[] batch_opm_a;
            delete[] batch_opm_b;
        }
        else
        { // for receiver-side memory
            delete[] batch_ot_rev_msg;
            delete[] batch_opm_rev_memory;
            delete[] batch_opm_delta;
        }
        delete[] batch_column_sum;
    }
};

#endif