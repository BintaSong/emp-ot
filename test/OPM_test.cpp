#include "emp-ot/emp-ot.h"
#include "test/test.h"
#include "test/OPM.h" 
#include <cstdlib>

using namespace std;

int port, party;
const static int threads = 2;

void test(int height, int index){

	// OPM<NetIO> rev(ALICE, threads, height);
	unique_ptr<OPM<NetIO>> rev_ptr( new  OPM<NetIO>(ALICE, threads, height));

	//rev.test_ggm_tree(height, index);
	// rev_ptr->test_opm();
	//rev_ptr->tes_opv(height);
	rev_ptr->test_opm();
	//std::cout<<"test() end."<<std::endl;

}

int main(int argc, char** argv){
	// int long long len = 2147483648;
	// std::cout<<"len = " << len << ", size = " << len*256/8/1000000 <<std::endl;
	// block *fuck = new block[len];
	// delete[] fuck;

	// block *fuck = new block[128];
	// block *pos = fuck + 100;
	// for(int i = 0; i < 29; i++) {
	// 	pos[i] = zero_block; 
	// }
	// delete fuck; 

	test(atoi(argv[1]), atoi(argv[2]));
	//test(5, 1);
	std::cout<<"main() end."<<std::endl;
}