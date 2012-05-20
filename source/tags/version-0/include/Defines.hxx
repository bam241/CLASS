#ifndef DEFINES_H
#define DEFINES_H 1

//#define DEBUG

#ifdef DEBUG
#define DBGL    cout << __FILE__ << ":" << __LINE__ << "[" << __FUNCTION__ << "]" << endl;
#define DBGV(x) cout << __FILE__ << ":" << __LINE__ << "[" << __FUNCTION__ << "]" #x  "='" << x << endl;
#else
#define DBGL
#define DBGV(x)
#endif


#endif
