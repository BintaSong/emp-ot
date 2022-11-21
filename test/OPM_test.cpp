#include "emp-ot/emp-ot.h"
#include "test/test.h"
#include "test/OPM.h" 
#include <cstdlib>

using namespace std;

int port, party;
const static int threads = 2;

void test_OPM(int height, int index){

	OPM<NetIO> rev(ALICE, threads, 8);
	// rev.receiver_init();
	// rev.random_permutation();
	// for(auto it = rev.permutation.begin() ; it != rev.permutation.end(); ++it) {
    //         std::cout<< *it << ' ';
	// }
	// std::cout<<std::endl; 

	rev.test_ggm_tree(height, index);

}

int main(int argc, char** argv){
	test_OPM(atoi(argv[1]), atoi(argv[1])); 
	// test_OPM(8, 4);  
}
