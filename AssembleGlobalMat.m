function [K, f] = AssembleGlobalMat(para, K, f, elem, elem_node, EK, Ef)
  for ii = 1 : para.nodes_per_elem
    m = elem_node(elem, ii);
    f(m) = f(m) + Ef(ii);
    for jj = 1 : para.nodes_per_elem
      n = elem_node(elem, jj);
      K(m, n) = K(m, n) + EK(ii, jj);
    end  % jj
  end  % ii
end
