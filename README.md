# Multi-agent optimal transport

This project implements a scalable, distributed algorithm for the optimal transport of large-scale multi-agent systems. 
Formulated as a problem of steering the collective towards a target probability measure while minimizing the total cost of transport
and implemented using a distributed online optimization algorithm to respect information constraints imposed by the multi-agent setting.  
Using optimal transport theory, the solution is realized as an iterative transport based on a stochastic proximal descent scheme. 
At each stage of the transport, the agents implement an online, distributed primal-dual algorithm to obtain local estimates of the Kantorovich potential for optimal transport from
the current distribution of the collective to the target distribution. 
Using these estimates as their local objective functions, the agents then implement the transport by stochastic proximal descent. 
This two-step process is carried out recursively by the agents to converge asymptotically to the target distribution. <br />

For more details, refer to: [preprint of paper](https://arxiv.org/pdf/1804.01572)

## Citation

If you found this work useful, please cite:
```@article{krishnan2025distributed,
  title={Distributed online optimization for multi-agent optimal transport},
  author={Krishnan, Vishaal and Mart{\'\i}nez, Sonia},
  journal={Automatica},
  volume={171},
  pages={111880},
  year={2025},
  publisher={Elsevier}
}```
