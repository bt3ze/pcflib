#ifndef OTEXTENSION_H_
#define OTEXTENSION_H_

#include "OTExtension/util/typedefs.h"
#include "OTExtension/util/crypto/crypto.h"
#include "OTExtension/util/socket.h"
//#include "OTExtension/ot/iknp-ot-ext-snd.h"
//#include "OTExtension/ot/iknp-ot-ext-rec.h"
#include "OTExtension/ot/alsz-ot-ext-snd.h"
#include "OTExtension/ot/alsz-ot-ext-rec.h"
//#include "../ot/nnob-ot-ext-snd.h"
//#include "../ot/nnob-ot-ext-rec.h"
#include "OTExtension/util/cbitvector.h"
#include "OTExtension/ot/xormasking.h"
#include "OTExtension/util/rcvthread.h"
#include "OTExtension/util/sndthread.h"
#include "OTExtension/util/channel.h"
#include "OTExtension/util/parse_options.h"

#include <vector>
#include <sys/time.h>

#include <limits.h>
#include <iomanip>
#include <string>

#include "Bytes.h"

using namespace std;

//TODO only for debugging purpose!!
static const char* m_cConstSeed[2] = {"437398417012387813714564100", "15657566154164561"};

USHORT		m_nPort = 7766;
const char* m_nAddr ;// = "localhost";

BOOL OT_Init(crypto* crypt);
BOOL OT_Cleanup();
BOOL OT_Connect();
BOOL OT_Listen();

void InitOTSender(const char* address, int port, crypto* crypt);
void InitOTReceiver(const char* address, int port, crypto* crypt);

BOOL ObliviouslyReceive(CBitVector& choices, CBitVector& ret, int numOTs, int bitlength, snd_ot_flavor stype, rec_ot_flavor rtype, crypto* crypt);
BOOL ObliviouslySend(CBitVector& X1, CBitVector& X2, int numOTs, int bitlength, snd_ot_flavor stype, rec_ot_flavor rtype, crypto* crypt);

// Network Communication
CSocket* m_vSocket;
uint32_t m_nPID; // thread id
field_type m_eFType;
uint32_t m_nBitLength;
MaskingFunction* m_fMaskFct;

// Naor-Pinkas OT
//BaseOT* bot;
OTExtSnd *sender;
OTExtRec *receiver;

SndThread* sndthread;
RcvThread* rcvthread;

uint32_t m_nNumOTThreads;
uint32_t m_nBaseOTs;
uint32_t m_nChecks;

bool m_bUseMinEntCorAssumption;
ot_ext_prot m_eProt;

double rndgentime;

int32_t read_test_options(int32_t* argcp, char*** argvp, uint32_t* role, uint64_t* numots, uint32_t* bitlen,
		uint32_t* secparam, string* address, uint16_t* port, ot_ext_prot* protocol, snd_ot_flavor* sndflav,
		rec_ot_flavor* rcvflav, uint32_t* nthreads, uint32_t* nbaseots, uint32_t* nchecks, bool* usemecr, uint32_t* runs);


void OT_alsz_send(const char* addr, unsigned short port, uint64_t num_OTs, uint32_t bitlength, uint32_t sec_param, std::vector<Bytes>& send_vals1, std::vector<Bytes> & send_vals2);
void OT_alsz_recv(const char* addr, unsigned short port, uint64_t num_OTs, uint32_t bitlength, uint32_t sec_param, Bytes selection_bits, std::vector<Bytes> &result_bytes);


#endif
