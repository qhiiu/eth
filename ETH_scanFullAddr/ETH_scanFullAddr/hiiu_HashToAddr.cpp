#include <stdint.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <iostream>
#include <vector>
#include <cstdint>

using namespace std;

//----------------------------------------------------------------------------------
std::string hash160keccakToAddr(uint32_t* hash)
{
	char tmp[3];
	std::string addr;

	addr.append("0x");
	for (int i = 0; i < 20; i++) {
		sprintf(tmp, "%02x", ((uint8_t*)hash)[i]);
		addr.append(tmp);
	}

	return addr;
}

