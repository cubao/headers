// should sync
// - https://github.com/cubao/polyline-ruler/blob/master/src/cubao_inline.hpp
// - https://github.com/cubao/headers/tree/main/include/cubao/cubao_inline.hpp

#ifdef CUBAO_INLINE
#undef CUBAO_INLINE
#endif

#ifndef CUBAO_STATIC_LIBRARY
#define CUBAO_INLINE inline
#else
#define CUBAO_INLINE
#endif
