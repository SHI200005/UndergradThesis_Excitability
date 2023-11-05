# Excitability in Gene Regulatory Networks

[Summary Slides](https://github.com/SHI200005/UndergradThesis_Excitability/blob/main/annex/Eng_copyright.pdf).

Copyright © [2023] [L Shi @ Nanjing University]. All rights reserved.

The text, images, and their arrangement on this [repository] are all subject to copyright and other intellectual property protection. These materials may not be copied for use without prior written permission from [L Shi @ Nanjing University].

This [repository] may also contain materials that are subject to the copyright rights of [L Shi, the supervisor and the PhD student @ Nanjing University].

For inquiries regarding the use of materials on this [repository], please contact [slnsinlangmc@sina.com].

By accessing or using this [repository], you agree to abide by the terms of this Copyright Notice.

Last updated: [2023/08/07]

[L Shi @ Nanjing University]

## Information

title = {基因调控网络中的可兴奋性研究}, title* = {Excitability in Gene Regulatory Networks},

author* = {L Shi},

keywords = {系统生物学,基因调控网络,可兴奋性,噪声}, keywords* = {{Systems Biology},{Gene Regulatory Networks},Excitability,Noise},

grade = {2018}, department = {物理学院}, department* = {School of Physics}, major = {物理学}, major* = {Physics},

supervisor = {  ,教授}, supervisor*= {Professor   }, (the author may need extra permission to show this)

submit-date = {2022-06-07},

## Abstract

可兴奋性是系统面对刺激快速进入兴奋状态并在刺激结束后回到初始状态的能力。可兴奋性广泛存在于生物系统中，但是缺少系统性的研究。我们通过建立最简单的包含耦合的正负反馈环的三节点网络，来系统性地研究实现可兴奋性的条件与其对应的动力学。通过构建网络结构对应的常微分方程组，模拟并分析网络的动力学，我们发现12种包含正负反馈的三节点网络结构中，正反馈涉及三个节点或者两个节点同时促进另一个节点的六个网络结构难以实现可兴奋，其余则能够实现可兴奋性。考虑到内外噪声对于生物系统的影响，我们发现在部分网络结构中，系统兴奋后会产生“不应期”，不对外噪声响应；还有网络结构受到内噪声的影响，其兴奋状态呈现出类似于“计时器”的多个振荡脉冲现象。我们筛选出了实现可兴奋性的网络结构，帮助人们进一步了解可兴奋性，对有可兴奋性基因网络的合成有指导意义。

Excitability is the ability of systems to enter an excited state quickly in response to a stimulus and to return to the initial state when the stimulus is withdrawn. Excitability is widespread in biological systems but has not been systematically studied. We systematically investigate the conditions under which excitability is achieved and their corresponding kinetics by building the minimal three-node network comprising coupled positive and negative feedback loops. By constructing a system of ordinary differential equations and simulating and analysing the dynamics of the network, we found that of 12 three-node network structures, six network structures with positive feedback involving three nodes or two nodes simultaneously promoting another node are difficult to achieve excitability, while the others can allow for excitability. Considering the effect of internal and external noise on biological systems, we found that in some network structures, excitability is followed by a refractory period where the system does not respond to external noise. Some internal noise affects other network structures, and their excitable state is characterised by multiple oscillatory pulses similar to a genetic timer. The current thesis identifies the conditions for excitable genetic circuits, helping to further our understanding of excitability and promoting the synthesis of excitable gene networks.

## Computation tools

[MATLAB](https://www.mathworks.com/products/matlab.html) and [Oscill8](https://oscill8.sourceforge.net/).

## 插图目录

1.1 含有正负反馈耦合回路的基因网络

1.2 可兴奋神经元

1.3 枯草芽孢杆菌感受态

1.4 可兴奋性的单参数分岔图表述

---

2.1 含有正负反馈耦合回路的基因网络

2.2 单参数分岔图与行为分类——非可兴奋

2.3 单参数稳态图与行为分类——可兴奋

---

3.1 Lotka-Volterra 竞争模型

3.2 兴奋与稳定流形

3.3 可兴奋网络与不可兴奋网络

3.4 可兴奋网络参数取值频率

3.5 参数敏感性

---

4.1 B.subtilis 感受态模型

4.2 参数波动的阈值

4.5 图4.2(e) 不应期

---

5.1 Gillespie 方法

5.2 合成生物学手段鉴别内外噪声

5.3 参数与内外噪声下的敏感性的一致性

5.4 转录-翻译速率与内噪声

5.5 转录-翻译相对速率与内噪声(1)

5.6 转录-翻译相对速率与内噪声(2)

5.7 内噪声下的基因计时器

5.8 基因计时器

---

A.1 网络#9

A.2 网络#9 结合反应
