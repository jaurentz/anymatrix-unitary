function alpha = szego_verblunsky(mu, count)
%SZEGO_VERBLUNSKY Compute interior Verblunsky coefficients of a measure.
%   ALPHA = SZEGO_VERBLUNSKY(MU,K) returns alpha_0,...,alpha_(K-1).
%   No finite unitary boundary phase is included; supply BETA separately.
%   Currently supports SZEGO_JACOBI_MEASURE descriptors.
alpha = js_verblunsky(mu,count);
end
