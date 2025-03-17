#ifndef __DEFINES__
#define __DEFINES__

enum Gate_type {
    H,
    T,
    CX,
    SWAP
};

enum Circuit_type {
    Main,
    If,
    While
};
static const char* Gate_view[4] = {"h ", "t ", "cx ", "swap "};
#endif