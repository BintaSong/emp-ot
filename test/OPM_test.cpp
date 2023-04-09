#include "emp-ot/emp-ot.h"
#include "test/test.h"
#include "test/OPM.h"
#include "test/shared_shuffle.h"
#include <cstdlib>

using namespace std;

int port, party;
const static int threads = 2;

// template <typename IO>
// void batch_opt_gen_and_store(int party, int height, int n_opt, IO *io, OT<IO> *ot)
// {

// 	std::string filename = party == ALICE ? OPT_SENDER_FILE : OPT_RECEIVER_FILE;

// 	std::ofstream outfile(filename);
// 	if (outfile.is_open())
// 		outfile.close();
// 	else
// 		error("create a directory to store opm data");
// 	FileIO fio(filename.c_str(), false);
// 	fio.send_data(&party, sizeof(int));
// 	fio.send_data(&n_opt, sizeof(int));

// 	for (int i = 0; i < n_opt; i++)
// 	{
// 		// std::cout<< "\nbatch gen, for i = " << i <<std::endl;
// 		OPM<NetIO> opm(party, io, ot, 2, height);
// 		opm.test_opm_with_OT();
// 		if (party == ALICE)
// 		{ // ALICE is the OPM sender
// 			fio.send_data(opm.opm_a, sizeof(block) * opm.dimension);
// 			fio.send_data(opm.opm_b, sizeof(block) * opm.dimension);
// 		}
// 		else
// 		{
// 			fio.send_data(opm.permutation.data(), sizeof(int) * opm.dimension);
// 			fio.send_data(opm.opm_delta, sizeof(block) * opm.dimension);
// 		}
// 	}
// 	std::cout << "batch_opt_gen_and_store() DONE." << std::endl;
// }

void read_opt_test(int dimension, int n_opt)
{
	int in_party;
	// std::cout<< "before open " <<std::endl;
	FileIO rec_fio(OPT_RECEIVER_FILE.c_str(), true);
	FileIO sen_fio(OPT_SENDER_FILE.c_str(), true);
	// std::cout<< "before read " <<std::endl;
	rec_fio.recv_data(&in_party, sizeof(int));
	rec_fio.recv_data(&n_opt, sizeof(int));
	sen_fio.recv_data(&in_party, sizeof(int));
	sen_fio.recv_data(&n_opt, sizeof(int));

	OPT sen_opt(ALICE, dimension), rev_opt(BOB, dimension);
	block *opt_check = new block[dimension];
	for (int i = 0; i < n_opt; i++)
	{
		rec_fio.recv_data(rev_opt.permutation.data(), sizeof(int) * dimension);
		rec_fio.recv_data(rev_opt.delta, sizeof(block) * dimension);

		sen_fio.recv_data(sen_opt.a, sizeof(block) * dimension);
		sen_fio.recv_data(sen_opt.b, sizeof(block) * dimension);

		// check OPT relation
		for (int i = 0; i < dimension; i++)
		{
			// delta = p(a)+b
			opt_check[i] = rev_opt.delta[i] ^ sen_opt.a[rev_opt.permutation[i]] ^ sen_opt.b[i];
			if (!cmpBlock(opt_check + i, &zero_block, 1))
			{
				log::error("opt check fails.\n");
			}
		}
	}
	delete[] opt_check;
}

int main(int argc, char **argv)
{

	int port, party, height, n_bits_width; // make sure all functions work for non-power-of-two lengths
	parse_party_and_port(argv, &party, &port);
	height = atoi(argv[3]);
	n_bits_width = atoi(argv[4]);

	bool pre_setup = false;
	pre_setup = atoi(argv[5]); 

	NetIO *io = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port);
	IKNP<NetIO> *iknp = new IKNP<NetIO>(io, true);

	size_t width, dimension;
	dimension = 1 << (height - 2);
	width = 1 << (n_bits_width);


	// auto t = clock_start();
	// offline_opt_batch_gen(party, dimension, width, io, iknp);
	// std::cout << "\n" << n_bits_width << " " << height << " time for offline_opt_batch_gen(): " << time_from(t)<< " microseconds"  << std::endl;

	if (!pre_setup) {
		log::info("before OPT generation...\n");
		auto t = clock_start();
		batch_opt_gen_and_store(party, height, n_bits_width, io, iknp);
		std::cout << "\ntime for batch_opt_gen_and_store(): " << time_from(t)<< " microseconds"  << std::endl;
		log::info("OPT generation done.\n");
	}


	//This step will initilize GBN-based shuffle
	// unique_ptr<SharedShuffle<NetIO>> shuffle_ptr(new SharedShuffle<NetIO>(party, io, iknp, dimension, width));

	// // perform GBN shuffle
	// shuffle_ptr->test_shuffle_GBN(height-2, n_bits_width);
	// log::info("GBN shuffle done.\n");

	// test io latency
// 	int m = 100,n = 10000, s = n/m;
// 	block *share = new block[n];
// 	block *rec = new block[n];
// 	memset(share, 0, sizeof(block)*n);
// 	memset(rec, 0, sizeof(block)*n);

// std::cout <<" counter = " << io->counter  << std::endl;

// auto t2 = clock_start();
// 	if(party == ALICE) {
// 		//for(int i = 0; i < 10000; i++)
// 		io->send_block(share, n);
// 		io->recv_block(rec, n);
// 	}
// 	else{
// 		//for(int i = 0; i < n; i++)
// 		io->recv_block(rec, n);
// 		io->send_block(share, n);
// 	}
// 	io->flush();
// std::cout <<" t2 = " << time_from(t2)<< " microseconds"  << std::endl;

// sleep(2);

// std::cout <<" counter = " << io->counter  << std::endl;


// auto t1 = clock_start();

// 	if(party == ALICE) {
		
// 		for(int i = 0; i < n; i+=m){
// 		//std::cout <<" test t1 ..."  << std::endl;
// 		io->send_block(share+i, m);
// 		io->recv_block(rec+i, m);
// 		}
// 	}
// 	else{
// 		for(int i = 0; i < n; i+=m){
// 		 //std::cout <<" test alice t1 ..."  << std::endl;
// 			io->recv_block(rec+i, m);
// 			io->send_block(share+i, m);
// 		}
// 	}
// 	io->flush();
// std::cout <<" t1 = " << time_from(t1)<< " microseconds"  << std::endl;


// std::cout <<" counter = " << io->counter  << std::endl;


// 	delete[] share;
// 	delete[] rec;


	delete io;
	delete iknp;
	std::cout << "main() end." << std::endl;
}