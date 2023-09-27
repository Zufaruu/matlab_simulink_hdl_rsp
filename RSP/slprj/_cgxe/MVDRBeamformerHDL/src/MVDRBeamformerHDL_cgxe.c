/* Include files */

#include "MVDRBeamformerHDL_cgxe.h"
#include "m_Kwnlh4DDKo0d1E3PVche9F.h"

unsigned int cgxe_MVDRBeamformerHDL_method_dispatcher(SimStruct* S, int_T method,
  void* data)
{
  if (ssGetChecksum0(S) == 48542131 &&
      ssGetChecksum1(S) == 56547143 &&
      ssGetChecksum2(S) == 390480318 &&
      ssGetChecksum3(S) == 965225324) {
    method_dispatcher_Kwnlh4DDKo0d1E3PVche9F(S, method, data);
    return 1;
  }

  return 0;
}
