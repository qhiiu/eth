#include <stdint.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <iostream>
#include <vector>
#include <cstdint>

using namespace std;

void decodeAddrToHash160keccak(const char* addr, uint32_t* hash160Keccak)
{
	string address = addr;

	std::vector<unsigned char> hashORxpoint;
	// address.erase(0, 2); // 
	for (int i = 0; i < 40; i += 2) {
		uint8_t c = 0;
		for (size_t j = 0; j < 2; j++) {
			uint32_t c0 = (uint32_t)address[i + j];
			uint8_t c2 = (uint8_t)strtol((char*)&c0, NULL, 16);
			if (j == 0)
				c2 = c2 << 4; 
			c |= c2;
		}
		hashORxpoint.push_back(c);
	}
	
	for (size_t i = 0; i < hashORxpoint.size(); i++) {
		((uint8_t*)hash160Keccak)[i] = hashORxpoint.at(i);
	}
}
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// int main()
// {
// 	uint32_t h[5];
// 	decodeAddrToHash160keccak("0x73B1D0D7EEa322bAc6cfE7B9329daf79ca1Ac76A",h);
// 	for (int i = 0; i < 5; i++){	printf("\n --- h[] : %d ", h[i]); }

// }