/*-******************************************************* -- HEADER -{{{1- */
/*-	Named template arguments.
 *
 * See Boost.Parameter or the implementation by Babtiste Wicht [1] for some
 * discussion around the concept.
 *
 * tl;dr: Allows named template arguments as follows:
 *
 *   using ActualClassA = SomeTemplate<>; // use defaults
 *   using ActualClassB = SomeTemplate<   // override one specific option
 *   	ScalarType<float>
 *   >;
 *   using ActualClassC = SomeTemplate<   // order doesn't matter.
 *   	ResultType<double>,
 *   	Method<EMethod::forward>
 *   >;
 *
 *
 * The setup for the above code could be as follows:
 *
 *   template< typename tType > 
 *   using ScalarType = named::TypeArgument< struct ScalarType_, double, tType >;
 *   template< typename tType > 
 *   using ResultType = named::TypeArgument< struct ResultType_, long, tType >;
 *   
 *   template< EMethod tValue >
 *   using Method = named::ValueArgument< struct ValueArg_, EMethod, EMethod::backward, tValue >;
 *
 *   template< class... tArgs >
 *   struct SomeTemplate
 *   {
 *   	using scalar_type = named::get_type_t<ScalarType, tArgs...>;
 *      ...
 *
 *      static constexpr EMethod method = named::get_value_v<EMethod, Method, tArgs...>;
 *      ...
 *   };
 *
 * (get_value_v is available when __cplusplus is >= 201402, otherwise use the old
 * style get_value<>::value.)
 *
 *
 * Note a few quirks:
 *   - unknown template arguments are ignored, so SomeTemplate< int > in the
 *     above will compile (and the `int` is ignored).
 *
 *   - arguments are matched on their tag, so in the following
 *     
 *       struct Tag {};
 *
 *       template< typename tType >
 *       using SomeArg = named::TypedArgument< Tag, double, tType >;
 *       template< typename tType >
 *       using OtherArg = named::TypedArgument< Tag, double, tType >;
 *
 *     SomeArg and OtherArg will behave identically (which should not come as a 
 *     surprise -- they are the same type in C++).
 * 
 *   - only the first argument of a certain kind (=with a specific tag) is
 *     used: SomeTemplate< ScalarType<int>, ScalarType<float> > will just see
 *     the first ScalarType<>, and thus have scalar_type == int.
 *
 * [1] https://baptiste-wicht.com/posts/2015/03/named-optional-template-parameters-compile-time.html
 */
/*-***************************************************************** -}}}1- */

#ifndef NAMED_HPP_98F4C1D9_EB24_46D8_822F_A7EB4C6FA9F0
#define NAMED_HPP_98F4C1D9_EB24_46D8_822F_A7EB4C6FA9F0

namespace named
{
	template< class tTag, typename tDefault, typename tType = tDefault >
	struct TypeArgument
	{
		using Tag = tTag;
		using Type = tType;
		using Def = tDefault;
	};

	template< class tTag, typename tType, tType tDefault, tType tValue = tDefault >
	struct ValueArgument
	{
		using Tag = tTag;
		using Type = tType;
		
		static constexpr Type def = tDefault;
		static constexpr Type value = tValue;
	};


	template< template<typename> class tArg, typename... tArgs >
	struct get_type
	{
		using Tag_ = typename tArg<void>::Tag;
		using Def_ = typename tArg<void>::Def;

		template< typename, typename... > struct Impl_;
		template< typename, typename, typename, typename... > struct Impl0_;

		template< typename tT0, typename tHead, typename... tTail >
		struct Impl_<tT0,tHead,tTail...>
			: Impl0_<tT0,typename tHead::Tag,tHead,tTail...>
		{};
		template< typename tT0 >
		struct Impl_<tT0>
		{
			using type = Def_;
		};

		template< typename tT0, typename, typename, typename... tTail >
		struct Impl0_
			: Impl_<tT0, tTail...>
		{};
		template< typename tTag, typename tHead, typename... tTail >
		struct Impl0_<tTag,tTag,tHead,tTail...>
		{ 
			using type = typename tHead::Type;
		};

		using type = typename Impl_<Tag_, tArgs...>::type;
	};

	template< template<typename> class tArg, typename... tArgs >
	using get_type_t = typename get_type<tArg,tArgs...>::type;


	template< typename tType, template<tType> class tArg, typename... tArgs >
	struct get_value
	{
		static constexpr tType zero_ = tType{};
		using Tag_ = typename tArg<zero_>::Tag;

		using Type_ = typename tArg<zero_>::Type; //TODO-static_assert is_same<>
		static constexpr Type_ def_ = tArg<zero_>::def;

		template< typename, typename... > struct Impl_;
		template< typename, typename, typename, typename... > struct Impl0_;

		template< typename tT0, typename tHead, typename... tTail >
		struct Impl_<tT0,tHead,tTail...>
			: Impl0_<tT0,typename tHead::Tag,tHead,tTail...>
		{};
		template< typename tT0 >
		struct Impl_<tT0>
		{
			static constexpr Type_ value = def_;
		};

		template< typename tT0, typename, typename, typename... tTail >
		struct Impl0_
			: Impl_<tT0, tTail...>
		{};
		template< typename tTag, typename tHead, typename... tTail >
		struct Impl0_<tTag,tTag,tHead,tTail...>
		{ 
			static constexpr Type_ value = tHead::value;
		};

		static constexpr Type_ value = Impl_<Tag_,tArgs...>::value;
	};

#	if __cplusplus >= 201402
	template< typename tType, template<tType> class tArg, typename... tArgs >
	constexpr typename tArg<tType{}>::Type get_value_v = get_value<tType,tArg,tArgs...>::value;
#	endif // ~ C++14
}

#endif // NAMED_HPP_98F4C1D9_EB24_46D8_822F_A7EB4C6FA9F0

