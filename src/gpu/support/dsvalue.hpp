/*-******************************************************* -- HEADER -{{{1- */
/*-	Static / Dynamic Value
 *
 * TODO: constructor/assignment with type conversion?
 */
/*-***************************************************************** -}}}1- */

#ifndef DSVALUE_HPP_583BB35C_21B9_4AA1_A46E_9F691C28AD7F
#define DSVALUE_HPP_583BB35C_21B9_4AA1_A46E_9F691C28AD7F

#include <cassert>
#include <type_traits>

/** Static (compile-time) value
 *
 * `StaticValue` models a value \a tValue of type \a tType that is known at
 * compile time. `StaticValue` is designed to be used together with
 * `DynamicValue` and facilitates (relatively) transparent propagation of
 * compile-time values when these are known, with fall-backs to run-time values
 * when the values are not known on before-hand.
 *
 * Example:
   \code
   template< int tValue >
   Result do_computation( StaticValue<int> const& aValue );
   Result do_computation( DynamicValue<int> const& aValue );
  
   template< class tCount >
   Result f( tCount aCount )
   {
   	return do_computation( aCount );
   }
   
  
   template< class tCount >
   class C
   {
   	...
   	public:
   		void g( tCount aCount ) { f( aCount ); }
   };
  
   ...
  
   using Cstatic = C< StaticValue<int,5> >;
   using Cdynamic = C< DynamicValue<int> >;
  
   CStatic cs;
   cs.g( 5 ); // OK. Asserts that the value is indeed 5
  
   Cdynamic cd;
   cd.g( 5 ); // OK. Accepts any int.
   \endcode
 *
 * A `StaticValue` may be implicitly constructed from a value of \a tType;
 * however, it is an error if the this value does not match the static value \a
 * tValue. Common arithmetic operations with only `StaticValue`s as arguments
 * will typically result in `StaticValue`s.
 *
 * \see DynamicValue
 */
template< typename tType, tType tValue >
struct StaticValue
{
	using type = tType;
	static constexpr bool isStatic = true;
	static constexpr tType staticValue = tValue;
	
#	if __cplusplus >= 201402
	constexpr
	StaticValue( tType aValue = tValue )
	{
		assert( aValue == tValue );
	}
#	else // !C++14
	constexpr
	StaticValue( tType = tValue )
	{}
#	endif // ~ C++14

	StaticValue( StaticValue&& ) = default;
	StaticValue& operator= (StaticValue&&) = default;

	StaticValue( StaticValue const& ) = default;
	StaticValue& operator= (StaticValue const&) = default;

	constexpr operator tType() const { return staticValue; }
	constexpr tType value() const { return staticValue; }
};

/** Dynamic (run-time) value
 *
 * `DynamicValue` models a value of \a tType that is only known at run-time; it
 * is a light-weight wrapper of underlying type \a tType. `DynamicValue` is
 * designed to be used together with `StaticValue` and facilitates (relatively)
 * transparent propagation of compile-time values when these are known, with
 * fallbacks to run-time values when the values are not known on before-hand.

 * `DynamicValue`s cannot be default constructed, but must be passed a value
 * of \a tType (but they are copy-/move-constructible).
 *
 * See `StaticValue` for a usage example.
 *
 * \see StaticValue
 */
template< typename tType >
struct DynamicValue
{
	using type = tType;
	static constexpr bool isStatic = false;

	constexpr
	DynamicValue( tType aValue )
		: mValue(aValue)
	{}

	DynamicValue( DynamicValue&& ) = default;
	DynamicValue& operator= (DynamicValue&&) = default;

	DynamicValue( DynamicValue const& ) = default;
	DynamicValue& operator= (DynamicValue const&) = default;

	constexpr operator tType() const { return mValue; }
	constexpr tType value() const { return mValue; }

	private:
		tType mValue;
};


// meta functions
template< class tValue >
struct is_dynamic_value
	: std::false_type
{};
template< typename tType >
struct is_dynamic_value< DynamicValue<tType> >
	: std::true_type
{};

template< class tValue >
struct is_static_value
	: std::false_type
{};
template< typename tType, tType tValue >
struct is_static_value< StaticValue<tType,tValue> >
	: std::true_type
{};


#if __cplusplus >= 201402
template< class tValue >
constexpr bool is_dynamic_value_v = is_dynamic_value<tValue>::vale;
template< class tValue >
constexpr bool is_static_value_v = is_static_value<tValue>::value;
#endif // ~ C++14


// arithm. operators
template< typename tType, tType tValue > constexpr
auto operator+ (StaticValue<tType,tValue> const&) -> StaticValue<tType,+tValue>
{
	return {};
}
template< typename tType, tType tValue > constexpr
auto operator- (StaticValue<tType,tValue> const&) -> StaticValue<tType,-tValue>
{
	return {};
}

template< typename tLTy, tLTy tLVal, typename tRTy, tRTy tRVal > constexpr
auto operator+ (StaticValue<tLTy,tLVal> const&, StaticValue<tRTy,tRVal> const&)
	-> StaticValue<decltype(tLVal+tRVal),tLVal+tRVal> { return {}; }
template< typename tLTy, tLTy tLVal, typename tRTy, tRTy tRVal > constexpr
auto operator- (StaticValue<tLTy,tLVal> const&, StaticValue<tRTy,tRVal> const&)
	-> StaticValue<decltype(tLVal-tRVal),tLVal-tRVal> { return {}; }
template< typename tLTy, tLTy tLVal, typename tRTy, tRTy tRVal > constexpr
auto operator* (StaticValue<tLTy,tLVal> const&, StaticValue<tRTy,tRVal> const&)
	-> StaticValue<decltype(tLVal*tRVal),tLVal*tRVal> { return {}; }
template< typename tLTy, tLTy tLVal, typename tRTy, tRTy tRVal > constexpr
auto operator/ (StaticValue<tLTy,tLVal> const&, StaticValue<tRTy,tRVal> const&)
	-> StaticValue<decltype(tLVal/tRVal),tLVal/tRVal> { return {}; }


template< typename tLTy, tLTy tLVal, typename tRTy, tRTy tRVal > constexpr
auto operator& (StaticValue<tLTy,tLVal> const&, StaticValue<tRTy,tRVal> const&)
	-> StaticValue<decltype(tLVal&tRVal),tLVal&tRVal> { return {}; }
template< typename tLTy, tLTy tLVal, typename tRTy, tRTy tRVal > constexpr
auto operator| (StaticValue<tLTy,tLVal> const&, StaticValue<tRTy,tRVal> const&)
	-> StaticValue<decltype(tLVal|tRVal),tLVal|tRVal> { return {}; }
template< typename tLTy, tLTy tLVal, typename tRTy, tRTy tRVal > constexpr
auto operator^ (StaticValue<tLTy,tLVal> const&, StaticValue<tRTy,tRVal> const&)
	-> StaticValue<decltype(tLVal^tRVal),tLVal^tRVal> { return {}; }
template< typename tLTy, tLTy tLVal, typename tRTy, tRTy tRVal > constexpr
auto operator<< (StaticValue<tLTy,tLVal> const&, StaticValue<tRTy,tRVal> const&)
	-> StaticValue<decltype(tLVal<<tRVal),(tLVal<<tRVal)> { return {}; }
template< typename tLTy, tLTy tLVal, typename tRTy, tRTy tRVal > constexpr
auto operator>> (StaticValue<tLTy,tLVal> const&, StaticValue<tRTy,tRVal> const&)
	-> StaticValue<decltype(tLVal>>tRVal),(tLVal>>tRVal)> { return {}; }


#endif // DSVALUE_HPP_583BB35C_21B9_4AA1_A46E_9F691C28AD7F

