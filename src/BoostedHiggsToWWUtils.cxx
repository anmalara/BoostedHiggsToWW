#include "UHH2/BoostedHiggsToWW/include/BoostedHiggsToWWUtils.h"

//================================================================================
bool XOR( bool a, bool b){
  return (!a&&b) || (!b && a);
}

//================================================================================

std::string itoa(int i)
{
    char res[10];
    sprintf(res, "%d", i);
    std::string ret(res);
    return ret;
}

