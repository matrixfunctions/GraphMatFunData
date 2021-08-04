function Z=expmpol_wrapper(A)
    % Disable norm estimation to make comparison more fair
    Z=expmpol(A,NaN,0);
end
