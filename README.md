# 计算动大作业

### Generalized_Alpha();

**Call**

- GetMass();获取集中质量阵

  - 针对梁单元编写，因此目前只能根据单元长度计算单元质量，平面单元后面再写。
  - 集中质量阵需要把质量安置在节点的每一个位移自由度上，因此需要设置单元自由度的排列方式，我用了$a^e=[x,y,z,\theta_x,\theta_y.\theta_z]$
  - 

- GetDamp(L);获取阻尼矩阵

  - L是工况编号

  - 阻尼矩阵有多种形式，可以看书P35，写成单元协调阻尼矩阵需要用梁单元的型函数，目前为了简单起见直接采用Rayleigh阻尼。

    

- 
