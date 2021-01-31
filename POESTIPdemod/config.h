#ifndef CONFIG_H
#define CONFIG_H

   #define USE_FLOATS 1 //0 for doubles
   
   #if USE_FLOATS==1
      #define DECIMAL_TYPE float
   #else
      #define DECIMAL_TYPE double
   #endif

   #if defined (__WIN32__)
      #include <windows.h>
      #ifndef ENABLE_VIRTUAL_TERMINAL_PROCESSING
         #define ENABLE_VIRTUAL_TERMINAL_PROCESSING 0x0004
      #endif
   #endif
#endif