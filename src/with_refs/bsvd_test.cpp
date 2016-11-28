#include "binmat.h"

int main(int argc, char** argv) {
  block_t a(~block_t(0));
  std::cout << bm_bitset(a << 5) << std::endl;
  std::cout << bm_bitset(a >> 3) << std::endl;
  return 0;
}
