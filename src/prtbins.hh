#ifndef PRT_BINS_HH
#define PRT_BINS_HH

template <typename... II,
          typename Bin, typename... Ax, typename Container, typename Filler>
std::enable_if_t<sizeof...(II)!=sizeof...(Ax)>
prtbins(std::ostream& o,
  ivanp::binner<Bin,std::tuple<Ax...>,Container,Filler>& h, II... ii
) {
  const auto& a = std::get<sizeof...(II)>(h.axes());
  for (unsigned i=1, n=a.nbins()+2; i<n; ++i) {
    for (unsigned s=0; s<sizeof...(II); ++s) o << "  ";
    o << "\033[35m[" << a.lower(i) << ',';
    if (i==n-1) o << "∞";
    else o << a.upper(i);
    o << ")\033[0m" << ( sizeof...(Ax)-sizeof...(II)!=1 ? '\n' : ' ' );
    prtbins(o,h,ii...,i);
  }
}

template <typename... II,
          typename Bin, typename... Ax, typename Container, typename Filler>
std::enable_if_t<sizeof...(II)==sizeof...(Ax)>
prtbins(std::ostream& o,
  ivanp::binner<Bin,std::tuple<Ax...>,Container,Filler>& h, II... ii
) { o << h.bin(ii...) << '\n'; }

template <typename... II,
          typename Bin, typename... Ax, typename Container, typename Filler>
std::ostream& operator<<(std::ostream& o,
  const ivanp::named_ptr<
    ivanp::binner<Bin,std::tuple<Ax...>,Container,Filler>
  >& h
) {
  o << "\033[32m" << h.name << "\033[0m\n";
  prtbins(o,*h);
  return o;
}

#endif
