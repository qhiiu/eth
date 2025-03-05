// // MAXX.SetBase10("115792089237316195423570985008687907852837564279074904382605163141518161494336"); // full
// MAXX.SetBase10("11579208923731619542357098500868790785283756427907490438260516314151"); // remove 10th-end

#include "Timer.h"
#include "KeyHunt.h" 
#include "Base58.h"
#include <string>
#include <cassert>
#include <signal.h>
#include <unistd.h>
#include <iostream>
#include <cstdint>
#include <fstream>
#include <map>
#include <cstring>
#include "hiiu_decodeAddrToHash160keccak.cpp"

using namespace std;
 
#define RANDOM 0
#define INPUT 1

//-----------------------------------------------------------------------
void CtrlHandler(int signum) {
	printf("\nBYE");
	printf("\nBYE\n\n");
	exit(signum);
}
//-----------------------------------------------------------------------
void check_file_exist(){
	FILE* file;
    file = fopen("$.txt", "r");
    if (file!=NULL)  {   
        printf("\n $.txt exists ======$$$$====== \n");
    	printf("\n.\n.\n.\n.\n.\n --- DONE !! take your fucking money !! ---- \n.\n.\n.\n.\n.\n.\n");

        exit(-1);   
    }
}
//-----------------------------------------------------------------------
uint64_t check_data(std::string priv)
{
    std::string fileName;
    fileName = "xDataScanned.txt";

    ifstream file(fileName);
    string line;
	uint64_t n = 0;

    if (file.is_open()) {
        while (getline(file, line)) {
            n += 1;
            if (priv == line){
              std::cout << std::endl <<"-------- Had in Database !!! exit() -------"<< std::endl << std::endl;
              exit(0);
            }
        }
        file.close();
    } else {   
        std::cout << "\n\nErr file !!! Don't have file : " << fileName << std::endl;
        exit(-1);
    }
	return n;
}
//-----------------------------------------------------------------------

void init_value(int mode, uint64_t xN,Int& privDec, Int& rangeStart, Int& rangeEnd)
{
    Int MINN, MAXX;
    MINN.SetBase10("0");
    // MAXX.SetBase16("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364140"); // full
    // MAXX.SetBase10("115792089237316195423570985008687907852837564279074904382605163141518161494336"); // full
    MAXX.SetBase10("115792089237316195423570985008687907852837564279074904382605163141518"); // remove 9th-end
    // std::cout << std::endl << "MAXX-10 : " << MAXX.GetBase10();
    // std::cout << std::endl << "MAXX-16 : " <<  MAXX.GetBase16();

    switch (mode){
        case RANDOM:
            privDec.Rand(&MAXX);
            // std::cout << std::endl << "privDec-10 : " << privDec.GetBase10();
            // std::cout << std::endl << "privDec-16 : " << privDec.GetBase16();
            break;

        case INPUT:
            std::cout << "\nRANGE__INPUT : "<< MINN.GetBase10() << " - " << MAXX.GetBase10();
            char* input_privDec = new char[64];

            while (true) 
            {
                std::cout << "\ninput__privDec : "; 
                cin.getline(input_privDec, 64);

                privDec.SetBase10(input_privDec);     

                if (privDec.IsGreaterOrEqual(&MINN) && privDec.IsLowerOrEqual(&MAXX)){ 
                    break; 
                }
            }

            // break;
        }

        // set _1B
        Int _1B, _xN1B;
        _1B.SetBase10("1000000000"); 
        // set xN1B = _1B then multiple to xN
        _xN1B = _1B;
        _xN1B.Mult(xN); 


        //set --- rangeStart
        rangeStart = privDec; 
        rangeStart.Mult(&_1B);
        // std::cout << std::endl << "rangeStart-10 :  " << rangeStart.GetBase10();
        // std::cout << std::endl << "rangeStart-16 :  " << rangeStart.GetBase16() << std::endl;

        //set --- rangEnd
        rangeEnd = rangeStart;
        rangeEnd.Add(&_xN1B);
        // std::cout << std::endl << "rangeEnd-10 :    " << rangeEnd.GetBase10();
        // std::cout << std::endl << "rangeEnd-16 :    " << rangeEnd.GetBase16() << std::endl;

        //print privDec info 
        uint64_t nChecked = 0;
        Int privDec_copy = privDec; 
        for (int i = 0; i < xN; i++){
            std::cout << "\nprivDec ======> " << privDec_copy.GetBase10(); //print 
            nChecked = check_data(privDec_copy.GetBase10()); // check priv 
            privDec_copy.AddOne(); // increase priv 
        } 
        std::cout << "\n\nnChecked : " << nChecked ;    
}

//-----------------------------------------------------------------------
void test(Int& rangeStart, Int& rangeEnd){

    rangeStart.SetBase16("35da1daf9584d5308b54b0d753d93f59c5b55fdea71d8cccd9941ab730000000");
    rangeEnd.SetBase16("35da1daf9584d5308b54b0d753d93f59c5b55fdea71d8cccd9941ab740000000");

}
//-----------------------------------------------------------------------

// Function to trim leading and trailing spaces and newlines
void trim(std::string &str) {
    str.erase(0, str.find_first_not_of(" \t\r\n"));  // Erase leading whitespaces and newline characters
    str.erase(str.find_last_not_of(" \t\r\n") + 1);  // Erase trailing whitespaces and newline characters
}

//-----------------------------------------------------------------------

// Funtion to convert listAddr to array of data h[5]
void run(){

    check_file_exist(); // check file $.txt

	// Global Init
	Timer::Init();
	rseed(Timer::getSeed32());

	// bool gpuAutoGrid = true;
	vector<int> gpuId = { 0 };
	vector<int> gridSize;

	Int privDec, rangeStart, rangeEnd;

    //set value 
    int xN = 1;
    int mode = RANDOM; 
    init_value(mode, xN, privDec, rangeStart, rangeEnd);

    // // test-here -----------------------------
    // test(rangeStart, rangeEnd); 
    // // test-here - end -----------------------

	std::string outputFile = "$.txt";
    std::cout << "\n\nOUTPUT FILE  : " << outputFile;


    // =========================== listAddr -> arrData ---- start ===========================
    // std::string name_file_data = "data/eth_100_.txt";
    std::string name_file_data = "data/eth_1000_.txt";
    // std::string name_file_data = "data/eth_1000000_.txt";
    // std::string name_file_data = "data/test_data.txt"; // test

    std::cout << "\nName_file_data : " << name_file_data;
    ifstream file_data(name_file_data); 

    string addr;

    // count number of addrETH
    uint32_t n_addrETH = 0;

    if (file_data.is_open()) {
        while (getline(file_data, addr)) {
            n_addrETH++;
        }   
        file_data.close();
    } else {
        std::cout << "\n\n !!!!!!!!! Err file_data !!!!!!!!! exit() \n\n" << std::endl;   
        exit(-1); 
    }

    // init list arrDataETH
    uint32_t arrDataETH[5 * n_addrETH + 1];
    arrDataETH[0] = n_addrETH;

	uint32_t gen_hash160keccak[5];

    // reset n_counter
    n_addrETH = 0;

    ifstream fileData(name_file_data); 
    string addrLine;

    // load data : each addr => hash160keccak => save to arrDataETH
    if (fileData.is_open()) {
        while (getline(fileData, addrLine)) {

            trim(addrLine);  // remove extra spaces or newlines

            decodeAddrToHash160keccak(addrLine.c_str(), gen_hash160keccak);
            arrDataETH[5 * n_addrETH + 1] = gen_hash160keccak[0];
            arrDataETH[5 * n_addrETH + 2] = gen_hash160keccak[1];
            arrDataETH[5 * n_addrETH + 3] = gen_hash160keccak[2];
            arrDataETH[5 * n_addrETH + 4] = gen_hash160keccak[3];
            arrDataETH[5 * n_addrETH + 5] = gen_hash160keccak[4];

            n_addrETH++;

        }   
        fileData.close();

    } else {   std::cout << "Err file_data !!!" << std::endl;   }

    // =========================== listAddr -> arrData ---- end ===========================

    if (gridSize.size() == 0) {
        for (int i = 0; i < gpuId.size(); i++) {
            gridSize.push_back(-1);
            gridSize.push_back(128);
        }
	}

	signal(SIGINT, CtrlHandler);

	KeyHunt* keyHunt;
    bool should_exit = false;
	keyHunt = new KeyHunt(arrDataETH, outputFile, rangeStart, rangeEnd, privDec, xN, should_exit);
	keyHunt->Search(gpuId, gridSize, should_exit);

	delete keyHunt;
}

int main(){
    run();
    
	return 0;
};
