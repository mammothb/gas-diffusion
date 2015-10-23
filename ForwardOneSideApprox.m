function value = ForwardOneSideApprox(h, u_now, u_next, u_next_next)
  value = (-3 * u_now + 4 * u_next - u_next_next) / 2 / h;
end
