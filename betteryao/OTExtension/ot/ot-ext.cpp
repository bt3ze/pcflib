/**
 \file 		ot-extension.cpp
 \author 	michael.zohner@ec-spride.de
 \copyright __________________
 \brief		Methods for the OT Extension routine
 */

#include "ot-ext.h"


#ifdef FIXED_KEY_AES_HASHING
void FixedKeyHashing(AES_KEY_CTX* aeskey, BYTE* outbuf, BYTE* inbuf, BYTE* tmpbuf, uint64_t id, uint32_t bytessecparam, crypto* crypt) {
#ifdef HIGH_SPEED_ROT_LT
	((uint64_t*) tmpbuf)[0] = id ^ ((uint64_t*) inbuf)[0];
	((uint64_t*) tmpbuf)[1] = ((uint64_t*) inbuf)[1];
#else
	memset(tmpbuf, 0, AES_BYTES);
	memcpy(tmpbuf, (BYTE*) (&id), sizeof(int));

	for (int i = 0; i < bytessecparam; i++) {
		tmpbuf[i] = tmpbuf[i] ^ inbuf[i];
	}
#endif

	crypt->encrypt(aeskey, outbuf, tmpbuf, AES_BYTES);

#ifdef HIGH_SPEED_ROT_LT
	((uint64_t*) outbuf)[0] ^= ((uint64_t*) inbuf)[0];
	((uint64_t*) outbuf)[1] ^= ((uint64_t*) inbuf)[1];
#else
	for (int i = 0; i < bytessecparam; i++) {
		outbuf[i] = outbuf[i] ^ inbuf[i];
	}
#endif
}
#endif



