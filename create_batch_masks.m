% chop up whole-brain mask into batches

clear all;

[mask, V] = ccnl_load_mask('masks/mask.nii');

all_inds = find(mask);

batch_size = 1000;

b = 0;
for s = 1:batch_size:length(all_inds)
    e = min(length(all_inds), s + batch_size - 1);
    inds = all_inds(s:e);
    b = b + 1;

    filename = sprintf('masks/mask_batchsize=%d_batch=%d.nii', batch_size, b);
    filename

    V.fname = filename;

    bmask = zeros(size(mask));
    bmask(inds) = 1;
    spm_write_vol(V, bmask);
end
