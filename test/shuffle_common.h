
#ifndef OPM_CONSTANT_H__
#define OPM_CONSTANT_H__
#include <iostream>

static std::string OPM_KEY_SENDER_FILE = "./data/opm_key_sender";
static std::string OPM_KEY_RECEIVER_FILE = "./data/opm_key_receiver";
static std::string OPM_SENDER_FILE = "./data/opm_sender";
static std::string OPM_RECEIVER_FILE = "./data/opm_receiver";


static std::string OPT_SENDER_FILE = "./data/opt_sender";
static std::string OPT_RECEIVER_FILE = "./data/opt_receiver";


static size_t MAX_BATCH_OPM_SIZE = 10000;


class log {
public:
    static void info(std::string s) {
        std::cout<< s;
    }
    static void green(std::string s) {
        std::cout<<"\033[1;32m" + s + "\033[0m";
    }
    static void error(std::string s){
        std::cout<<"\033[1;31mERROR: " + s + "\033[0m";
    }
};
#endif