/*-******************************************************* -- HEADER -{{{1- */
/*-	2D Matrix
 *
 * TODO: maybe use a fixed size array when possible (SizeType = StaticValue...)?
 */
/*-***************************************************************** -}}}1- */

#ifndef MATRIX_HPP_192A3248_B449_499F_A0DC_AEE81D3B30A9
#define MATRIX_HPP_192A3248_B449_499F_A0DC_AEE81D3B30A9

#include <iterator>
#include <type_traits>

#include <cstddef>

#include "compat.hpp"
#include "dsvalue.hpp"
#include "dsvector.hpp"

namespace aspect
{
	struct MatrixRowMajor;
	struct MatrixColMajor;
}

namespace detail
{
	template< typename tIter, typename tConstIter >
	class MatElementView final
	{
		static_assert( std::is_same<typename std::iterator_traits<tIter>::iterator_category,std::random_access_iterator_tag>::value, "MatElementView: requries random access iterator." );
		
		public:
			constexpr MatElementView( tIter, tIter );

		public:
			using iterator = tIter;
			using const_iterator = tConstIter;

			using value_type = typename std::iterator_traits<tIter>::value_type;
			
		public:
			CONSTEXPR_14 value_type& operator[] (std::size_t);
			constexpr value_type const& operator[] (std::size_t) const;
		
		public:
			constexpr bool empty() const;
			constexpr std::size_t size() const;

			CONSTEXPR_14 tIter begin(), end();
			constexpr tConstIter begin() const, end() const;

		private:
			tIter mBeg, mEnd;
	};
}

template< typename tScalar, class tLayout = aspect::MatrixRowMajor,  class tNSizeType = DynamicValue<std::size_t>, class tMSizeType = DynamicValue<std::size_t> >
class Matrix final
{
	// note: not a literal type (see ~Matrix()), so the "constexpr" on the
	// member functions is probably not really that helpful?
	
	public:
		using value_type = tScalar;

		using NSizeType = tNSizeType;
		using MSizeType = tMSizeType;

		using Layout = tLayout;

		using ColIterator = typename Layout::template ColIterator<tScalar*,NSizeType,MSizeType>;
		using RowIterator = typename Layout::template RowIterator<tScalar*,NSizeType,MSizeType>;
		using ConstColIterator = typename Layout::template ColIterator<tScalar const*,NSizeType,MSizeType>;
		using ConstRowIterator = typename Layout::template RowIterator<tScalar const*,NSizeType,MSizeType>;

		using ColView = detail::MatElementView<ColIterator,ConstColIterator>;
		using RowView = detail::MatElementView<RowIterator,ConstRowIterator>;
		using ConstColView = detail::MatElementView<ConstColIterator,ConstColIterator>;
		using ConstRowView = detail::MatElementView<ConstRowIterator,ConstRowIterator>;

		using ElementView = detail::MatElementView<tScalar*,tScalar const*>;
		using ConstElementView = detail::MatElementView<tScalar const*,tScalar const*>;

		using SizeType = decltype(std::declval<NSizeType>()*std::declval<MSizeType>());
		
	public:
		constexpr Matrix();
		constexpr Matrix( tScalar*, NSizeType, MSizeType );

		Matrix( NSizeType, MSizeType );
		Matrix( NSizeType, MSizeType, tScalar );

		~Matrix();

		Matrix( Matrix const& );
		Matrix& operator= (Matrix const&);

		Matrix( Matrix&& ) = default;
		Matrix& operator= (Matrix&&) = default;


	public:
		CONSTEXPR_14 tScalar& operator[] (std::size_t);
		constexpr tScalar const& operator[] (std::size_t) const;

		CONSTEXPR_14 tScalar& operator() (std::size_t,std::size_t);
		constexpr tScalar const& operator() (std::size_t,std::size_t) const;

	public:
		CONSTEXPR_14 tScalar* data();
		constexpr tScalar const* data() const;

		constexpr NSizeType n() const;
		constexpr MSizeType m() const;
		constexpr SizeType size() const;

		CONSTEXPR_14 RowView row(std::size_t);
		constexpr ConstRowView row(std::size_t) const;

		CONSTEXPR_14 ColView col(std::size_t);
		constexpr ConstColView col(std::size_t) const;

		CONSTEXPR_14 ElementView linear_elements();
		constexpr ConstElementView linear_elements() const;

		CONSTEXPR_14 tScalar* lend();
		CONSTEXPR_14 tScalar* lbegin();

		constexpr tScalar const* lend() const;
		constexpr tScalar const* lbegin() const;
		
	private:
		tScalar* mData;
		
		NSizeType mN;
		MSizeType mM;
};

#include "detail/matrix.inl"
#endif // MATRIX_HPP_192A3248_B449_499F_A0DC_AEE81D3B30A9
