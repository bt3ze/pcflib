#ifndef YAO_H_
#define YAO_H_

#include "YaoBase.h"
#include "garbled_circuit.h"

class Yao : public YaoBase
{

public:
	Yao(EnvParams &params);
	virtual ~Yao() { }

	virtual void start();

private:
	//void oblivious_transfer();
	void circuit_evaluate();
	void proc_gen_out();
	void proc_evl_out();

	// variables for Yao protocol
        std::vector<Bytes>          m_gen_inp_masks;
        std::vector<Bytes>          m_rnds;
	//vector<GarbledCct>     m_ccts;
        std::vector<garbled_circuit_t> m_gcs;

};

#endif
