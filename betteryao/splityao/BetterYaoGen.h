#ifndef BETTERYAOGEN_H_
#define BETTERYAOGEN_H_

#include "../YaoBase.h"
#include "BetterYao5.h"
#include "../garbled_circuit_m.h"

#include <vector>

class BetterYaoGen : public BetterYao5
{
public:
	BetterYaoGen(EnvParams &params);
	virtual ~BetterYaoGen() {}

	void start();

	//void oblivious_transfer();
	//void cut_and_choose();
	//void cut_and_choose2();
	void consistency_check();
	void circuit_evaluate();

private:
	void ot_init();
	void ot_random(); // sender has m pairs of l-bit strings, and receiver has m bits
	void cut_and_choose2_ot();
	void cut_and_choose2_precomputation();
	void cut_and_choose2_evl_circuit(size_t ix);
	void cut_and_choose2_chk_circuit(size_t ix);

	Bytes flip_coins(size_t len);

	void proc_gen_out();
	void proc_evl_out();

};

#endif /* BETTERYAOGEN_H_ */
