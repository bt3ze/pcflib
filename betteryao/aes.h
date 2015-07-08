// aes.h
// just a couple of function prototypes
// see aes.ccp for implementation

#ifndef AES_H_
#define AES_H_

#include "Bytes.h"

void KDF128_Fixed_Key(const uint8_t *in, uint8_t * out);

Bytes KDF128(const Bytes &in, const Bytes &key);
Bytes KDF256(const Bytes &in, const Bytes &key);

void KDF128(const uint8_t *in, uint8_t *out, const uint8_t *key);
void KDF256(const uint8_t *in, uint8_t *out, const uint8_t *key);

#ifdef __CPLUSPLUS
extern "C"
{
#endif

void AES_128_Key_Expansion(const uint8_t *userkey, uint8_t *key_schedule);
void AES_192_Key_Expansion(const uint8_t *userkey, uint8_t *key_schedule);
void AES_256_Key_Expansion(const uint8_t *userkey, uint8_t *key_schedule);

void AES_ECB_decrypt (const uint8_t *in, uint8_t *out, unsigned long length, const uint8_t *KS, int nr);
void AES_ECB_encrypt (const uint8_t *in, uint8_t *out, unsigned long length, const uint8_t *KS, int nr);

#ifdef __CPLUSPLUS
};
#endif



#endif
