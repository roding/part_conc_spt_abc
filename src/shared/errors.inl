namespace error
{
	template< typename... tArgs > inline
	RuntimeError::RuntimeError( char const* aFmt, tArgs&&... aArgs )
		: std::runtime_error( tfm::format( aFmt, std::forward<tArgs>(aArgs)... ) )
	{}
}
