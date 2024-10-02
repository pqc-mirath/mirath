//
// Finite field arithmetic, including vector and matrix operations
//
#ifndef MIRATH_ARITH_H
#define MIRATH_ARITH_H


#if (defined(_Ib_) || defined(_IIIb_) || defined(_Vb_))
#include "q_2/data_type.h"
#include "q_2/ff.h"
#include "q_2/matrix_ff_arith.h"
#include "q_2/matrix_ff_mu.h"
#include "q_2/vector_ff_arith.h"
#include "q_2/vector_ff_mu.h"

#elif (defined(_Ia_) || defined(_IIIa_) || defined(_Va_))
#include "q_16/data_type.h"
#include "q_16/ff.h"
#include "q_16/mu_fast/ff_mu.h"
//#include "q_16/mu_short/ff_mu.h"
#include "q_16/matrix_ff_arith.h"
#include "q_16/matrix_ff_mu.h"
#include "q_16/vector_ff_arith.h"
#include "q_16/vector_ff_mu.h"

//#elif (defined(_I?_) || defined(_III?_) || defined(_V?_))
//

#else
// By default we take q = 16
#include "q_16/data_type.h"
#include "q_16/ff.h"
#include "q_16/matrix_ff_arith.h"
#include "q_16/matrix_ff_mu.h"
#include "q_16/vector_ff_arith.h"
#include "q_16/vector_ff_mu.h"
#endif


#endif // MIRATH_ARITH_H