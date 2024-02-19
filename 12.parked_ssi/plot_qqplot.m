x0 = [-3:.1:3];
y = normpdf(x0,0,1);
plot(x0,y);hold on
histogram(xs,100,'Normalization','pdf')
qqplot(xs)


pwelch(squeeze(dataset_test(:,1,[1,2,3,13,14,15])),126,[],[],20)
legend('1','2','3','4','5','6')