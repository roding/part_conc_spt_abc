/*-******************************************************* -- HEADER -{{{1- */
/*-	C++ version compatibility
 */
/*-***************************************************************** -}}}1- */

#ifndef COMPAT_HPP_5EB8B4C2_58EF_49F4_83A6_D21D359A2DCF
#define COMPAT_HPP_5EB8B4C2_58EF_49F4_83A6_D21D359A2DCF

#include <memory>

//-/
//-/ C++14 language - relaxed constexpr
#if __cplusplus >= 201402
#	define CONSTEXPR_14 constexpr
#else
#	define CONSTEXPR_14 /*nothing*/
#endif

//-/
//-/ C++17 language - [[deprecated]] attribute
#if __cplusplus >= 201703
#	define ATTR_DEPRECATED [[deprecated]]
#else
#	if defined(__GNUC__) && __GNUC__ >= 6
#		define ATTR_DEPRECATED __attribute__((deprecated))
#	endif

#	if !defined(ATTR_DEPRECATED)
#		define ATTR_DEPRECATED /*nothing*/
#	endif
#endif


//-/
//-/ C++14 stdlib - make_unique()
#if __cplusplus >= 201402
using std::make_unique;
#else
template< typename tType, typename... tArgs > inline
std::unique_ptr<tType> make_unique( tArgs&&... aArgs )
{
	return std::unique_ptr<tType>( new tType( std::forward<tArgs>(aArgs)... ) );
}
#endif

//-/
//-/ C++14 stdlib - selected *_t variable templates
#if __cplusplus >= 201402
using std::enable_if_t;
#else
template< bool B, class T = void >
using enable_if_t = typename std::enable_if<B,T>::type;
#endif

#endif // COMPAT_HPP_5EB8B4C2_58EF_49F4_83A6_D21D359A2DCF
