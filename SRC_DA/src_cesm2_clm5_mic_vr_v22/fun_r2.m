function r2 = fun_r2(obs, mod)
    r2 = 1 - sum((mod - obs).^2, 'omitnan')/sum((obs - mean(obs, 'omitnan')).^2, 'omitnan');
end
