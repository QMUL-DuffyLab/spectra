using DelimitedFiles
using LinearAlgebra
using Formatting
using PyPlot

function CD_calc(dir)
  mu_file   = format("{}{}", dir, "mu_site.out")
  com_file  = format("{}{}", dir, "c_o_m.dat")
  eig_file  = format("{}{}", dir, "eigvecs.out")
  en_file   = format("{}{}", dir, "eigvals.out")
  mu        = readdlm(mu_file)
  r         = readdlm(com_file)
  eig       = readdlm(eig_file)
  en        = readdlm(en_file)
  ns = 2048
  nchl = 14
  λ = 1.0
  CD = zeros(ns)
  wn = zeros(ns)
  g_sum = zeros(ns)
  chi = Array{Real}(undef, nchl, ns, 2)
  for k = 1:nchl
    chi_file    = format("{}{}{:02d}{}", dir, "chi_i_", k, ".dat")
    chi[k, :, :]      = readdlm(chi_file) 
  end
  wn = chi[1, :, 1]

  for n = 1:nchl
    for m = 1:nchl
      dd = ((mu[n, :]) × (mu[m, :])) ⋅ (r[n, :] - r[m, :])
      g_sum = zeros(ns)
      for k = 1:nchl
        # chi_file = format("{}{}{:02d}{}", dir, "chi_i_", k, ".dat")
        # chi      = readdlm(chi_file) 
        # if (n == 1 && m == 1 && k == 1)
        #   wn       = chi[:, 1]
        # end
        # s        = sum(chi[:, 2])Cross product - W
        g        = chi[k, :, 2] / sum(chi[k, :, 2])
        g_sum .+= (g .* (eig[k, n] * eig[k, m] * en[k]))
      end
      CD .+= (g_sum .* dd)
    end
  end

  CD .*= (-16 * pi^2 / 9λ) .* wn
  # CD .*= (-16 * pi^2 / 9λ) .* wn
  wl = replace!(1e7 ./ wn, Inf=>NaN)

  open(format("{}{}", dir, "CD.dat"), "w") do io
    writedlm(io, [wl CD])
  end

  return [wl wn CD]
end

dir = ARGS[1]
if (dir[end] != "/")
  dir *= "/"
end
res = CD_calc(dir)

# p = plot(res[:, 1], res[:, 3], xaxis = ("frequency (nm)", (580, 720)))
p = plot(res[:, 1], res[:, 3])
ax = gca()
ax.set_xlim([580, 720])
xlabel("frequency (nm)")
ylabel("CD")
savefig(format("{}{}", dir, "CD.pdf"))

