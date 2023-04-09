#ifndef SHARED_SHUFFLE_H__
#define SHARED_SHUFFLE_H__

#include "OPM.h"
static const int MAX_N_BITS_GBN_WIDTH = 64;

typedef struct OPT
{
    size_t party, dimension;

    std::vector<size_t> permutation;
    block *delta = nullptr, *a = nullptr, *b = nullptr;

    OPT(int party, size_t dimension)
    {
        this->party = party;
        this->dimension = dimension;

        if (party == ALICE)
        {
            this->a = new block[dimension];
            this->b = new block[dimension];
        }
        else
        {
            this->permutation.resize(dimension);
            this->delta = new block[dimension];
        }
    }
    // TODO:
    OPT(const OPT &opt)
    {
        this->party = opt.party;
        this->dimension = opt.dimension;

        if (party == ALICE)
        {
            this->a = new block[opt.dimension];
            memcpy(this->a, opt.a, sizeof(block) * dimension);
            this->b = new block[opt.dimension];
            memcpy(this->b, opt.b, sizeof(block) * dimension);
        }
        else
        {
            this->permutation.resize(opt.dimension);
            memcpy(this->permutation.data(), opt.permutation.data(), sizeof(size_t) * dimension);
            this->delta = new block[opt.dimension];
            memcpy(this->delta, opt.delta, sizeof(block) * dimension);
        }
    }

    ~OPT()
    {
        // log::info("~OPT()\n");
        if (delta != nullptr)
            delete[] delta;
        if (a != nullptr)
            delete[] a;
        if (b != nullptr)
            delete[] b;
    }
} OPT;

template <typename IO>
void batch_opt_gen_and_store(int party, int height, int n_bits_width, IO *io, OT<IO> *ot)
{
    size_t width, dimension;
    width = 1 << n_bits_width;
    dimension = 1 << (height - 2);
    int n_opt = ceil(2.0 * n_bits_width / (height - 2) - 1) * (width / dimension);

    std::string filename = party == ALICE ? OPT_SENDER_FILE : OPT_RECEIVER_FILE;
    filename = filename + "_" + std::to_string(height - 2) + "_" + std::to_string(n_bits_width);

    std::ofstream outfile(filename);
    if (outfile.is_open())
        outfile.close();
    else
        error("create a directory to store opm data");
    FileIO fio(filename.c_str(), false);
    fio.send_data(&party, sizeof(int));
    fio.send_data(&n_opt, sizeof(int));

    for (int i = 0; i < n_opt; i++)
    {
        // std::cout<< "\nbatch gen, for i = " << i <<std::endl;
        OPM<NetIO> opm(party, io, ot, 2, height);
        opm.test_opm_with_OT();
        if (party == ALICE)
        { // ALICE is the OPM sender
            fio.send_data(opm.opm_a, sizeof(block) * opm.dimension);
            fio.send_data(opm.opm_b, sizeof(block) * opm.dimension);
        }
        else
        {
            fio.send_data(opm.permutation.data(), sizeof(size_t) * opm.dimension);
            fio.send_data(opm.opm_delta, sizeof(block) * opm.dimension);
        }
    }
    std::cout << "batch_opt_gen_and_store() DONE." << std::endl;
}

template <typename IO>
void offline_opt_batch_gen(int party, size_t dimension, size_t width, IO *io, OT<IO> *ot)
{
    std::cout << "offline_opt_batch_gen() begin..." << std::endl;
    int n_bits_dimension, n_bits_width;
    n_bits_dimension = log2(dimension);
    n_bits_width = log2(width);

    int d = ceil(1.0 * n_bits_width / n_bits_dimension);
    int n_opt = (2 * d - 1) * (width / dimension);

    std::string filename = party == ALICE ? OPT_SENDER_FILE : OPT_RECEIVER_FILE;
    filename = filename + "_OFF_TEST_" + std::to_string(n_bits_dimension) + "_" + std::to_string(n_bits_width);

    std::ofstream outfile(filename);
    if (outfile.is_open())
        outfile.close();
    else
        error("create a directory to store opm data");
    FileIO fio(filename.c_str(), false);
    fio.send_data(&party, sizeof(int));
    fio.send_data(&n_opt, sizeof(int));

    std::cout << "offline_opt_batch_gen() will generate " << n_opt << " OPTs." << std::endl;

    for (int i = 0; i < n_opt; i++)
    {
        // std::cout<< "\nbatch gen, for i = " << i <<std::endl;
        OPM<NetIO> opm(party, io, ot, 2, n_bits_dimension + 2);
        opm.test_opm_with_OT();
        if (party == ALICE)
        { // ALICE is the OPM sender
            fio.send_data(opm.opm_a, sizeof(block) * opm.dimension);
            fio.send_data(opm.opm_b, sizeof(block) * opm.dimension);
        }
        else
        {
            fio.send_data(opm.permutation.data(), sizeof(size_t)  * opm.dimension);
            fio.send_data(opm.opm_delta, sizeof(block) * opm.dimension);
        }
    }
    std::cout << party << " : offline_opt_batch_gen() DONE." << std::endl;
}

template <typename IO>
class SharedShuffle
{
public:
    int party, threads_n;
    int n_bits_dimension, n_opt, current_opt = 0;

    // configuration for generalized benes network
    size_t dimension, gbn_width, n_bits_gbn_width, n_gbn_layer, n_pems_per_layer, gbn_d;

    std::string opt_filename;
    block *msg_buf = nullptr;

    IO *io;
    OT<IO> *ot; // NOTE: may not use

    /*
    Why use shared_pt here:
        OPT maintains some constructor-allocated heap storage, thus we define copy constructor
        because it will be used when opt->pushback() is called. Using shared_ptr can delay
        deconstruction of OPT to the memory deletion of opt.
    */
    std::vector<std::shared_ptr<OPT>> opt;
    // std::vector<OPT *> opt;

    void init_gbn(size_t dimension, size_t width)
    {
        this->dimension = dimension;
        this->gbn_width = width;

        this->n_bits_gbn_width = int(log2(width));
        if (1 << n_bits_gbn_width != width)
            log::error("gbn_width must be a power-of-two.");

        // if (int(log2(gbn_width)) % int(log2(dimension)) != 0)
        //     log::error("log2(dimension) does not divide log2(gbn_width)\n");
        this->n_gbn_layer = ceil(1.0 * log2(gbn_width) / log2(dimension));
        this->n_pems_per_layer = size_t(gbn_width / dimension);
        this->gbn_d = 2 * n_gbn_layer - 1;
    }

    SharedShuffle(int party, IO *io, OT<IO> *ot, size_t dimension, size_t gbn_width, int threads_n = 1, bool pre_setup = true)
    {
        this->party = party;
        this->io = io;
        this->ot = ot;
        this->dimension = dimension;
        this->n_bits_dimension = int(log2(dimension));
        if (1 << n_bits_dimension != dimension)
            log::error("opt dimension must be a power-of-two.");
        this->threads_n = threads_n;
        this->current_opt = 0;

        init_gbn(dimension, gbn_width);

        // TODO: assume pre_setup is done.
        opt_filename = (party == ALICE ? OPT_SENDER_FILE : OPT_RECEIVER_FILE);
        opt_filename = opt_filename + "_OFF_MALI_TEST_" + std::to_string(n_bits_dimension) + "_" + std::to_string(n_bits_gbn_width);
        read_opt_from_file(opt_filename);
        // if opts are not enough, then generate
        // if (gbn_d * n_pems_per_layer > n_opt) {
        //     int n_bits_dimension = int(log2(dimension));
        //     log::info("before OPT generation...\n");
        //     batch_opt_gen_and_store(party, n_bits_dimension+2, gbn_d * n_pems_per_layer - n_opt, io, ot);
        //     log::info("OPT generation done.\n");
        // }

        std::cout << "Shuffle requires " << gbn_d * n_pems_per_layer << " OPTs. " << n_opt << " are provided." << std::endl;
        if (n_opt < gbn_d * n_pems_per_layer)
        {
            // std::cout<< "Shuffle requires " << gbn_d * n_pems_per_layer << " OPTs. But only"<< n_opt << " are provided." <<std::endl;
            log::error("SharedShuffle(): No sufficient OPTs\n");
        }

        msg_buf = new block[n_opt*dimension];
    }

    ~SharedShuffle()
    {
        if (msg_buf != nullptr)
            delete[] msg_buf;
    }

    void read_opt_from_file(std::string filename)
    {
        FileIO fio(filename.c_str(), true);
        int in_party, n_opt;

        // NOTE: OPT file structure: [party][n_opt][n_opt OPTs]
        fio.recv_data(&in_party, sizeof(int));
        fio.recv_data(&n_opt, sizeof(int));
        this->n_opt = n_opt;
        if (in_party != party)
            error("wrong party");

        for (int i = 0; i < n_opt; i++)
        {
            std::shared_ptr<OPT> tmp(new OPT(party, dimension));
            // std::cout<< "for i = " << i <<std::endl;
            if (party == BOB)
            {
                // std::cout<< "BOB aaa : " << i <<std::endl;
                fio.recv_data(tmp->permutation.data(), sizeof(size_t) * dimension);
                // std::cout<< "BOB : " << i <<std::endl;
                fio.recv_data(tmp->delta, sizeof(block) * dimension);
            }
            else
            {
                // std::cout<< "Alice aaa : " << i <<std::endl;
                fio.recv_data(tmp->a, sizeof(block) * dimension);
                // std::cout<< "Alice : " << i <<std::endl;
                fio.recv_data(tmp->b, sizeof(block) * dimension);
            }
            // std::cout<< "before push " << i <<std::endl;
            opt.push_back(tmp);
        }
    }

    void sender_shuffle(block *share, block *out_msg)
    {
        if (current_opt >= n_opt)
        {
            log::error("sender_shuffle(): no OPTs for use\n");
            return;
        }
        for (size_t i = 0; i < dimension; i++)
        {
            // log::info("sender_shuffle "+std::to_string(i)+":\n");
            out_msg[i] = share[i] ^ opt[current_opt]->a[i];
            share[i] = opt[current_opt]->b[i];
        }
        current_opt++;
    }

    void receiver_shuffle(block *share, block *in_msg)
    {
        if (current_opt >= n_opt)
        {
            log::error("receiver_shuffle(): no OPTs for use\n");
            return;
        }
        for (size_t i = 0; i < dimension; i++)
        {
            // s[i] = s[p(i)] ^ in_msg[p(i)] ^ Delta[i]
            block *tmp = new block[dimension];
            memcpy(tmp, share, sizeof(block) * dimension);
            share[i] = tmp[opt[current_opt]->permutation[i]] ^ in_msg[opt[current_opt]->permutation[i]] ^ opt[current_opt]->delta[i];
        }
        current_opt++;
    }

    /*gbn_position_map() maps `opt_pos` to the corresponding position in GBN, using `layer` and `pem_index`.
        Do read the following help info:
            1. `layer` of gbn is indexed following [n_layer-1, n_layer-2,...,2,1,0,1,2,..., n_layer-2, n_layer-1]
            2. `std::bitset` encodes integer using [higher bits <-- lower bits] format.
                1) we want organize bits_ret_pos = [higher part of bit_gbn_pos] || bit_opt_pos || [lower part of bit_gbn_pso].
            3.
    */
    // size_t gbn_position_map(int layer, size_t pem_index, int opt_pos)
    // {
    //     /*
    //     Do read the following help info:
    //         1. `layer` of gbn is indexed following [n_layer-1, n_layer-2,...,2,1,0,1,2,..., n_layer-2, n_layer-1]
    //         2. `std::bitset` encodes integer using [higher bits <-- lower bits] format.
    //             1) we want organize bits_ret_pos = [higher part of bit_gbn_pos] || bit_opt_pos || [lower part of bit_gbn_pso].
    //         3.
    //     */
    //     std::bitset<MAX_N_BITS_GBN_WIDTH> bits_gbn_pos(pem_index), bits_opt_pos(opt_pos);

    //     // letf-shift `n_bit_dimension` bits for layer `i` in [`n_gbn_layer`-1, `layer`]
    //     for (int i = n_gbn_layer - 1; i >= layer; i--)
    //     {
    //         int start = i * n_bits_dimension;
    //         for (int j = 0; j < n_bits_dimension; j++)
    //         {
    //             bits_gbn_pos[start + j + n_bits_dimension] = bits_gbn_pos[start + j]; // copy bits from i-th layer to (i+1)-th
    //         }
    //     }
    //     // set bits in `layer` to be `bits_opt_pos`
    //     for (int j = 0; j < n_bits_dimension; j++)
    //     {
    //         bits_gbn_pos.set(layer * n_bits_dimension + j, bits_opt_pos[j]);
    //     }
    //     size_t ret = bits_gbn_pos.to_ulong();
    //     return ret;
    // }

    size_t gbn_position_map(int layer, size_t pem_index, int opt_pos)
    {
        /*
        Do read the following help info:
            1. `layer` of gbn is indexed following [n_layer-1, n_layer-2,...,2,1,0,1,2,..., n_layer-2, n_layer-1]
            2. `std::bitset` encodes integer using [higher bits <-- lower bits] format.
                1) we want organize bits_ret_pos = [higher part of bit_gbn_pos] || bit_opt_pos || [lower part of bit_gbn_pso].
            3.
        */
        std::bitset<MAX_N_BITS_GBN_WIDTH> bits_gbn_pos(pem_index), bits_opt_pos(opt_pos);

        // compute the [left, right] regin to be left-shifted
        int left, right;
        left = n_bits_gbn_width - n_bits_dimension; // len(`pem_indx`) = n_bits_gbn_width - n_bits_dimension
        right =  layer == 0 ? 0: n_bits_gbn_width - (n_gbn_layer-1)*n_bits_dimension; 

        // perform left-shift by `n_bits_dimension` bits
        for(int j = left; j < right; j++) {
            bits_gbn_pos[j+n_bits_dimension] = bits_gbn_pos[j];
        }
        
        // set `bits_opt_pos`
        for (int j = 0; j < n_bits_dimension; j++)
        {
            bits_gbn_pos.set(right + j, bits_opt_pos[j]);
        }

        size_t ret = bits_gbn_pos.to_ulong();
        return ret;
    }

    void sender_shuffle_gbn(block *share, block *out_msg, int layer_index, int pem_index)
    {
        if (current_opt >= n_opt)
        {
            log::error("sender_shuffle_gbn(): no OPTs for use\n");
            return;
        }
        for (size_t i = 0; i < dimension; i++)
        {
            size_t gbn_pos = gbn_position_map(layer_index, pem_index, i);
            out_msg[i] = share[gbn_pos] ^ opt[current_opt]->a[i];
            share[gbn_pos] = opt[current_opt]->b[i];
        }
        current_opt++;
    }

    void receiver_shuffle_gbn(block *share, block *in_msg, int layer_index, int pem_index)
    {
        if (current_opt >= n_opt)
        {
            log::error("receiver_shuffle_gbn(): no OPTs for use\n");
            return;
        }
        block *tmp = new block[dimension];
        for (size_t i = 0; i < dimension; i++)
        {
            size_t gbn_pos_source = gbn_position_map(layer_index, pem_index, opt[current_opt]->permutation[i]);
            tmp[i] = share[gbn_pos_source];
        }
        for (size_t i = 0; i < dimension; i++)
        {
            // s[i] = s[p(i)] ^ in_msg[p(i)] ^ Delta[i]
            size_t gbn_pos_dest = gbn_position_map(layer_index, pem_index, i);
            share[gbn_pos_dest] = tmp[opt[current_opt]->permutation[i]] ^ in_msg[opt[current_opt]->permutation[i]]; //^ opt[current_opt]->delta[i];
        }
        current_opt++;
        delete[] tmp;
    }

    void shuffle_with_network(block *share, block *msg)
    {
        if (party == ALICE)
        {
            sender_shuffle_gbn(share, msg);
            io->send_data(msg, dimension * sizeof(block));
            io->flush();
        }
        else
        {
            io->recv_data(msg, dimension * sizeof(block));
            receiver_shuffle_gbn(share, msg);
        }
    }
    /*
    `share` is with length `width`;
    `msg` is with length `opt_dimension`.
    */
    void shuffle_with_network_gbn(block *share, block *msg, int layer_index, int pem_index)
    {
        if (party == ALICE)
        {
            sender_shuffle_gbn(share, msg, layer_index, pem_index);
            io->send_data(msg, dimension * sizeof(block));
            io->flush();
        }
        else
        {
            io->recv_data(msg, dimension * sizeof(block));
            receiver_shuffle_gbn(share, msg, layer_index, pem_index);
        }
    }

    void shuffle_with_buf_gbn(block *share, block *msg, int layer_index, int pem_index)
    {
        if (party == ALICE)
        {
            sender_shuffle_gbn(share, msg, layer_index, pem_index);
            memcpy(&msg_buf[(current_opt - 1) * dimension], msg, dimension);
        }
        else
        {
            memcpy(msg, &msg_buf[current_opt * dimension], dimension);
            receiver_shuffle_gbn(share, msg, layer_index, pem_index);
        }
    }

    void middle_shuffle_with_buf_gbn(block *share, block *msg) //FIXME: msg is not used in this version 
    {
        if(party == ALICE)
        {
            for(size_t i = 0; i < n_pems_per_layer; i++) {
                for(size_t j = 0; j < dimension; j++) {
                    msg_buf[current_opt * dimension + j] = share[i*dimension+j] ^ opt[current_opt]->a[j];
                    share[i*dimension+j] = opt[current_opt]->b[j];
                }
                current_opt++;
            }
        }
        else {
            block *tmp = new block[dimension]; 
            for(size_t i = 0; i < n_pems_per_layer; i++) {
                memcpy(tmp, &share[i*dimension], dimension);
                for(size_t j = 0; j < dimension; j++){
                    share[i*dimension+j] = tmp[opt[current_opt]->permutation[j]] ^ msg_buf[i*dimension + opt[current_opt]->permutation[j]];
                }
                current_opt++;
            }
        }
    }

    // FIXME: so far it onshuffle_with_networkly works for special-case GBN: `n_bits_gbn_width` is dividable by `n_bits_dimension`!
    void simple_shuffle_with_network_gbn(block *share, block *msg)
    {
        // shuffle from top to down [n_gbn_layer-1 --> 1]
        for (int i = n_gbn_layer - 1; i > 0; i--)
        {
            for (int j = 0; j < n_pems_per_layer; j++)
            {
                shuffle_with_network_gbn(share, msg, i, j);
            }
        }

        for (int j = 0; j < n_pems_per_layer; j++)
        {
            shuffle_with_network_gbn(share, msg, 0, j);
        }
        //middle_shuffle_with_buf_gbn(share, msg);

        // shuffle from down to top [1 --> n_gbn_layer-1]
        for (int i = 1; i < n_gbn_layer - 1; i++)
        {
            for (int j = 0; j < n_pems_per_layer; j++)
            {
                shuffle_with_network_gbn(share, msg, i, j);
            }
        }
    }

    // FIXME: so far it onshuffle_with_networkly works for special-case GBN: `n_bits_gbn_width` is dividable by `n_bits_dimension`!
    // void simple_shuffle_with_buf_gbn(block *share, block *msg)
    // {
    //     // for BOB, receive gbn messages first
    //     if (party == BOB)
    //     {   
    //         // std::cout<< "I am BOB" <<std::endl;
    //         io->recv_block(msg_buf, n_opt * dimension);
    //         //io->send_block(msg_buf, 1);
    //     }

    //     // shuffle from top to down [n_gbn_layer-1 --> 1]
    //     for (int i = n_gbn_layer - 1; i > 0; i--)
    //     {
    //         for (int j = 0; j < n_pems_per_layer; j++)
    //         {
    //             shuffle_with_buf_gbn(share, msg, i, j);
    //         }
    //     }

    //     // shuffle the middle layer = 0
    //     // for (int j = 0; j < n_pems_per_layer; j++)
    //     // {
    //     //     shuffle_with_buf_gbn(share, msg, 0, j);
    //     // }
    //     middle_shuffle_with_buf_gbn(share, msg);

    //     // shuffle from down to top [1 --> n_gbn_layer-1]
    //     for (int i = 1; i < n_gbn_layer - 1; i++)
    //     {
    //         for (int j = 0; j < n_pems_per_layer; j++)
    //         {
    //             shuffle_with_buf_gbn(share, msg, i, j);
    //         }
    //     }

    //     // for ALICE, send gbn messages at last
    //     if (party == ALICE)
    //     {
    //         io->send_block(msg_buf, n_opt * dimension);
    //         io->flush();
    //         //io->send_block(msg_buf, 1);
    //     }
    // }

    void simple_shuffle_with_buf_gbn(block *share, block *msg)
    {
        // for BOB, receive gbn messages first
        if (party == BOB)
        {   
            // std::cout<< "I am BOB" <<std::endl;
            auto t = clock_start();
            io->recv_block(msg_buf, n_opt * dimension);
            
            //io->send_block(msg_buf, 1);
        }

        // shuffle from top to down [n_gbn_layer-1 --> 1]
        for (int i = n_gbn_layer - 1; i > 0; i--)
        {
            for (int j = 0; j < n_pems_per_layer; j++)
            {
                shuffle_with_buf_gbn(share, msg, i, j);
            }
        }

        // shuffle the middle layer = 0
        // for (int j = 0; j < n_pems_per_layer; j++)
        // {
        //     shuffle_with_buf_gbn(share, msg, 0, j);
        // }
        middle_shuffle_with_buf_gbn(share, msg);

        // shuffle from down to top [1 --> n_gbn_layer-1]
        for (int i = 1; i < n_gbn_layer - 1; i++)
        {
            for (int j = 0; j < n_pems_per_layer; j++)
            {
                shuffle_with_buf_gbn(share, msg, i, j);
            }
        }

        // for ALICE, send gbn messages at last
        if (party == ALICE)
        {
            io->send_block(msg_buf, n_opt * dimension);
            io->flush();
            //io->send_block(msg_buf, 1);
        }
    }

    /*
        void test_shuffle() {

            // first read opt data from local file
            std::string filename = (party == ALICE ? OPT_SENDER_FILE : OPT_RECEIVER_FILE);
            read_opt_from_file(filename);
            // std::cout<< "Read done." <<std::endl;

            // the prepared fake share and simulate local shuffle
            block *sen_share = new block[dimension], *rev_share = new block[dimension];
            memset(sen_share, 0, sizeof(block)*dimension);
            memset(rev_share, 0, sizeof(block)*dimension);

            // perform shuffle with network
            block *msg = new block[dimension];
            if(party == ALICE) {
                sender_shuffle(sen_share, msg);
                io->send_data(msg, dimension*sizeof(block)); io->flush();
            }
            else {
                io->recv_data(msg, dimension*sizeof(block));
                receiver_shuffle(rev_share, msg);
            }

            block *share = party == ALICE ? sen_share : rev_share;
            for(int i = 0; i < dimension; i++) {
                std::cout<< share[i] << " " <<std::endl;
            }
            std::cout<<std::endl;

            delete[] sen_share;
            delete[] rev_share;
            delete[] msg;
        }
    */

    void test_shuffle()
    {
        // using fake share only for test
        block *share = new block[dimension];
        memset(share, 0, sizeof(block) * dimension);

        // perform shuffle with network
        block *msg = new block[dimension];
        shuffle_with_network(share, msg);
        for (size_t i = 0; i < dimension; i++)
        {
            std::cout << share[i] << std::endl;
        }
        std::cout << std::endl;

        delete[] share;
        delete[] msg;
    }

    void test_shuffle_GBN(int n_bits_dimension, int n_bits_width)
    {
        // std::cout<< "Read done." <<std::endl;
        // only for test
        block *share = new block[gbn_width];
        memset(share, 1, sizeof(block) * gbn_width);

        // perform shuffle with GBN
        block *msg = new block[dimension];

        auto t2 = clock_start();
        simple_shuffle_with_network_gbn(share, msg);
        std::cout << "\ntime for simple_shuffle_with_network_gbn(): " << time_from(t2) << " microseconds" << std::endl;

        // for(int i = 0; i < gbn_width; i++){
        //     std::cout<< share[i] <<std::endl;
        // }
        // std::cout<<std::endl;

        delete[] share;
        delete[] msg;
    }

    void test_shuffle_GBN_with_buf(int n_bits_dimension, int n_bits_width)
    {
        // std::cout<< "Read done." <<std::endl;
        // only for test
        block *share = new block[gbn_width];
        memset(share, 1, sizeof(block) * gbn_width);

        // perform shuffle with GBN
        block *msg = new block[dimension];

        auto t2 = clock_start();
        simple_shuffle_with_buf_gbn(share, msg);
        std::cout << "n_bits_dimension: " << n_bits_dimension << ", n_bits_width: "<< n_bits_width << ". Time for simple_shuffle_with_buf_gbn(): " << time_from(t2) << " microseconds" << std::endl;

        delete[] share;
        delete[] msg;
    }

    // //TODO: compute bucket size from statistical security parameter and opm dimension
    // int get_bucket_size(int stat_sec = 40) {

    // }

    // //TODO: combine multiple OPTs as one single OPT
    // void opt_composition(opm *new_opt, const int size) {

    // }

    // //TODO: run shuffle using `opt + benes shuffle network` approach
    // void perform_shuffle_with_benes(opm *opt, const block *in, block *out) {

    // }

    // //TODO: run consistency check after shuffle
    // bool shuffle_consistency_check(const block *share) {

    // }

    // //TODO: SPDZ fake offline line. This is only used for test
    // void fake_spdz_offline() {

    // }
};

#endif
