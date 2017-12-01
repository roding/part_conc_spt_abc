/*-******************************************************* -- HEADER -{{{1- */
/*-	basic_instrument.hpp 
 *
 * Very basic helpers for very basic instrumentation/"profiling" via wall-time.
 */
/*-***************************************************************** -}}}1- */

#ifndef BASIC_INSTRUMENT_HPP_1D8771F9_665B_48A7_8D27_AA511A1FBF9D
#define BASIC_INSTRUMENT_HPP_1D8771F9_665B_48A7_8D27_AA511A1FBF9D

#include <chrono>
#include <cstdint>

using InstrTime = std::uint64_t;
using InstrClock = std::chrono::high_resolution_clock;

#define INSTR_BLOCK_BEGIN( ident ) auto const _instrBeg##ident = InstrClock::now()
#define INSTR_BLOCK_END( ident, accum ) do { \
		auto const _instrEnd##ident = InstrClock::now(); \
		accum += std::chrono::duration_cast<std::chrono::duration<InstrTime,std::micro>>(_instrEnd##ident - _instrBeg##ident).count(); \
	} while(0) \
	/*ENDM*/


#if INSTR_CURRENT_LEVEL >= 1
#	define INSTR_LEVEL1(...) __VA_ARGS__
#	define INSTR_BLOCK1_BEGIN(ident) INSTR_BLOCK_BEGIN(ident)
#	define INSTR_BLOCK1_END(ident,accum) INSTR_BLOCK_END(ident, accum)
#else
#	define INSTR_LEVEL1(...) (void)0
#	define INSTR_BLOCK1_BEGIN(ident) (void)0
#	define INSTR_BLOCK1_END(ident,accum) do {} while(0)
#endif // ~ INSTR_CURRENT_LEVEL

#if INSTR_CURRENT_LEVEL >= 2
#	define INSTR_LEVEL2(...) __VA_ARGS__
#	define INSTR_BLOCK2_BEGIN(ident) INSTR_BLOCK_BEGIN(ident)
#	define INSTR_BLOCK2_END(ident,accum) INSTR_BLOCK_END(ident, accum)
#else
#	define INSTR_LEVEL2(...) (void)0
#	define INSTR_BLOCK2_BEGIN(ident) (void)0
#	define INSTR_BLOCK2_END(ident,accum) do {} while(0)
#endif // ~ INSTR_CURRENT_LEVEL

#if INSTR_CURRENT_LEVEL >= 3
#	define INSTR_LEVEL3(...) __VA_ARGS__
#	define INSTR_BLOCK3_BEGIN(ident) INSTR_BLOCK_BEGIN(ident)
#	define INSTR_BLOCK3_END(ident,accum) INSTR_BLOCK_END(ident, accum)
#else
#	define INSTR_LEVEL3(...) (void)0
#	define INSTR_BLOCK3_BEGIN(ident) (void)0
#	define INSTR_BLOCK3_END(ident,accum) do {} while(0)
#endif // ~ INSTR_CURRENT_LEVEL

#if INSTR_CURRENT_LEVEL >= 4
#	define INSTR_LEVEL4(...) __VA_ARGS__
#	define INSTR_BLOCK4_BEGIN(ident) INSTR_BLOCK_BEGIN(ident)
#	define INSTR_BLOCK4_END(ident,accum) INSTR_BLOCK_END(ident, accum)
#else
#	define INSTR_LEVEL4(...) (void)0
#	define INSTR_BLOCK4_BEGIN(ident) (void)0
#	define INSTR_BLOCK4_END(ident,accum) do {} while(0)
#endif // ~ INSTR_CURRENT_LEVEL

#endif // BASIC_INSTRUMENT_HPP_1D8771F9_665B_48A7_8D27_AA511A1FBF9D
