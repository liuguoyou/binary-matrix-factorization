#include <cstdio>

int main() {
  unsigned char a = ~(unsigned char)0;
  unsigned char b = (a >> 1);
  unsigned long mask = 0x80;
  unsigned int i;
  for (i = 0; i < sizeof(unsigned char)*8; i++, mask >>=1) {
    putchar(mask & (a >> 1) ? '1': '0');
  }
  putchar('\n');
  return 0;
}
