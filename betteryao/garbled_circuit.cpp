#include "garbled_circuit.h"

Bytes setBytes( __m128i & i)
{

        Bytes q;
        q.resize(16,0);
        _mm_storeu_si128(reinterpret_cast<__m128i*>(&q[0]), i);

        return q;
}

void *copy_key(void *old_key)
{
	__m128i *new_key = 0;
	if (old_key != 0)
	{
		new_key = (__m128i*)_mm_malloc(sizeof(__m128i), sizeof(__m128i));
		*new_key = *reinterpret_cast<__m128i*>(old_key);
	}
	return new_key;
}

void delete_key(void *key)
{
	if (key != 0) _mm_free(key);
}
