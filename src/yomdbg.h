USE YOMDBG, ONLY : PCRC
#define crc(x) call pcrc(__FILE__, __LINE__, x, #x)
#define CRC(x) call pcrc(__FILE__, __LINE__, x, #x)
