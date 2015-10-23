function value = CentralApprox(h, u_next, u_prev)
  value = (u_next - u_prev) / 2 / h;
end
