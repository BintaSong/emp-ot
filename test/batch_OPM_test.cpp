#include "emp-ot/emp-ot.h"
#include "test/test.h"
#include "test/BatchOPM.h"
#include "test/shared_shuffle.h"
#include <cstdlib>

using namespace std;

int port, party;
const static int threads = 2;

int main(int argc, char **argv)
{

	int port, party, depth, n_bits_width; // make sure all functions work for non-power-of-two lengths
	parse_party_and_port(argv, &party, &port);
	depth = atoi(argv[3]);
	n_bits_width = atoi(argv[4]);

	bool pre_setup = false;
	pre_setup = atoi(argv[5]);

	bool malicious = false;
#ifdef MALICIOUS
	malicious = true;
#endif

	NetIO *io = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port);
	IKNP<NetIO> *iknp = new IKNP<NetIO>(io, malicious);

	size_t width, dimension, n_threads;
	dimension = 1 << depth;
	width = 1 << (n_bits_width);
	n_threads = 8; // atoi(argv[6]);

	if (!pre_setup)
	{
		BatchOPM<NetIO> bopm(party, io, iknp, n_threads, depth);
		size_t batch_size = bopm.get_gbn_batch_size(dimension, width);
#ifdef MALICIOUS
		bopm.test_batch_opm_with_OT(width, dimension, batch_size);
#else
		bopm.test_batch_semi_opm_with_OT(width, dimension, batch_size);
#endif
	}
	// This step will initilize GBN-based shuffle
	// if (n_bits_width % (height - 2) == 0)
	{
		unique_ptr<SharedShuffle<NetIO>> shuffle_ptr(new SharedShuffle<NetIO>(party, io, iknp, dimension, width));
		// perform GBN shuffle
		std::cout << "Begin online shuffle..." << std::endl;
		// shuffle_ptr->test_shuffle_GBN_with_buf(height - 2, n_bits_width);
		shuffle_ptr->test_shuffle_GBN_with_buf(depth, n_bits_width);
		log::info("Online shuffle done.\n");
	}
	// sleep(5);
	delete io;
	delete iknp;
	std::cout << "main() end." << std::endl;
}