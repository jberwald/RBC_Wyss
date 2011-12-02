% This program is used to generate a new movie with 1 ninth of original
% pixels. 
% written by Jia-Rong Yeh  on January 2009.

clear;
    %filnam=sprintf('C_4.mat',i);  % load the original movie
    filnam = sprintf('../New_3/new_3.txt')
    load(filnam);
    [new_bund,new_img]=cond_img(img_bund,img);  % generate the new movie
    filnam=sprintf('../New_3/new_3_cond.mat');
    save(filnam,'new_bund','new_img');       % save the new movie and boundary map
    clear img_bund img new_bund new_img;
 