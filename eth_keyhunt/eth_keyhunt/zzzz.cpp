
Int k(&key);
k.Add((uint64_t)incr);
Point p = secp->ComputePublicKey(&k);

std::string addr = secp->GetAddressETH(p); 

std::string addr = secp->GetAddressETH(hash);

