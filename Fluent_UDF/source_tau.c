#include "udf.h"
#include "prop.h"
#include "pdf_transport.h"
#include "flamelet.h"
#include "math.h"

/*source term for non reactive, non diffusive residence time*/
DEFINE_SOURCE(tau_source,c,t,dS,eqn)/*attach to the UDS-1 source term in fluid cell zone*/
{
  real source;
  dS[eqn] = 0;
  return 1.0;
}
