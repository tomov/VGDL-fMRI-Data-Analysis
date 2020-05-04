% little script to print code to output
% used for reproducibility (e.g. to know which version of the code produced which output)

S = dbstack();
p = S(2).file;
disp(['--------------------- BEGIN FILE ', p, '-------------------']);
code = fileread(p);
disp(code);
disp(['--------------------- END FILE ', p, '-------------------']);


