// aes.hs
// note that some of this file is borrowed from the JustGarble project, which borrows from openssl.
// this library links to openssl anyway, but some of the declarations here are different
// so we include both (TODO: clean this up)

// the following copyright from openssl is included

/* crypto/aes/aes.h -*- mode:C; c-file-style: "eay" -*- */
/* ====================================================================
 * Copyright (c) 1998-2002 The OpenSSL Project.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * 3. All advertising materials mentioning features or use of this
 *    software must display the following acknowledgment:
 *    "This product includes software developed by the OpenSSL Project
 *    for use in the OpenSSL Toolkit. (http://www.openssl.org/)"
 *
 * 4. The names "OpenSSL Toolkit" and "OpenSSL Project" must not be used to
 *    endorse or promote products derived from this software without
 *    prior written permission. For written permission, please contact
 *    openssl-core@openssl.org.
 *
 * 5. Products derived from this software may not be called "OpenSSL"
 *    nor may "OpenSSL" appear in their names without prior written
 *    permission of the OpenSSL Project.
 *
 * 6. Redistributions of any form whatsoever must retain the following
 *    acknowledgment:
 *    "This product includes software developed by the OpenSSL Project
 *    for use in the OpenSSL Toolkit (http://www.openssl.org/)"
 *
 * THIS SOFTWARE IS PROVIDED BY THE OpenSSL PROJECT ``AS IS'' AND ANY
 * EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE OpenSSL PROJECT OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 * ====================================================================
 *
 */


#ifndef AES_H_
#define AES_H_

#include <xmmintrin.h>              /* SSE instructions and _mm_malloc */
#include <emmintrin.h>              /* SSE2 instructions               */

#include "Bytes.h"


  // new fixed-key AES encryption
typedef __m128i block;
    
typedef struct { __m128i rd_key[15]; int rounds; } AES_KEY_J;

//void KDF128_Fixed_Key(const uint8_t *in, uint8_t * out);
//AES
void KDF128(uint8_t *out, const uint8_t * in, const AES_KEY_J * key);

// obsolete
//Bytes KDF128(const Bytes &in, const Bytes &key);
//Bytes KDF256(const Bytes &in, const Bytes &key);

//void KDF128(const uint8_t *in, uint8_t *out, const uint8_t *key);

//SHA
void KDF256(const uint8_t *in, uint8_t *out, const uint8_t *key);

#ifdef __CPLUSPLUS
extern "C"
{
#endif

  // old KSS AES encryption code
  // these functions were not actually supported
  // no one ever turned the AESNI flag on (with aes support)
//void AES_128_Key_Expansion(const uint8_t *userkey, uint8_t *key_schedule);
  //void AES_192_Key_Expansion(const uint8_t *userkey, uint8_t *key_schedule);
  //void AES_256_Key_Expansion(const uint8_t *userkey, uint8_t *key_schedule);

  //void AES_ECB_decrypt (const uint8_t *in, uint8_t *out, unsigned long length, const uint8_t *KS, int nr);
  //void AES_ECB_encrypt (const uint8_t *in, uint8_t *out, unsigned long length, const uint8_t *KS, int nr);


void AES_128_fixed_Key_Expansion(const uint8_t *userkey, uint8_t *key_schedule);
void AES_192_fixed_Key_Expansion(const uint8_t *userkey, uint8_t *key_schedule);
void AES_256_fixed_Key_Expansion(const uint8_t *userkey, uint8_t *key_schedule);


int AES_set_encrypt_key(const unsigned char *userKey, const int bits, AES_KEY_J *key);
void AES_set_decrypt_key_fast(AES_KEY_J *dkey, const AES_KEY_J *ekey);
int AES_set_decrypt_key(const unsigned char *userKey, const int bits, AES_KEY_J *key);

void AES_encrypt(const unsigned char *in, unsigned char *out, const AES_KEY_J *key);
void AES_decrypt(const unsigned char *in, unsigned char *out, const AES_KEY_J *key);
void AES_ecb_encrypt_blks(block *blks, unsigned nblks, AES_KEY_J *key);
void AES_ecb_decrypt_blks(block *blks, unsigned nblks, AES_KEY_J *key);

#ifdef __CPLUSPLUS
};
#endif



#endif
