---
layout: post
title: "AI Engineering Mindsets and Practices for Computational Biologists (Part I)"
date: 2026-05-25
description: "Exploring how Spec-Driven Development shifts AI-assisted coding from vibe coding to controlled system building, and why this engineering mindset matters deeply for computational biologists."
tags:
  - tech
  - engineering
  - sdd
categories:
  - tech-blog
bilingual: true
default_lang: en
related_posts: false
title_zh: "计算生物学家的AI工程思维和实践（一）"
description_zh: "探讨Spec-Driven Development如何将AI辅助编程从抽卡式写代码转变为可控的系统构建，以及这种工程思维为何对计算生物学家如此重要。"
---

<div
  class="bilingual-post"
  data-default-lang="en"
  data-title-en="AI Engineering Mindsets and Practices for Computational Biologists (Part I)"
  data-title-zh="计算生物学家的AI工程思维和实践（一）">
  <div class="bilingual-post-switch" role="group" aria-label="Language switcher">
    <button class="bilingual-post-button is-active" type="button" data-bilingual-target="en" aria-controls="ai-engineering-en" aria-pressed="true">English</button>
    <button class="bilingual-post-button" type="button" data-bilingual-target="zh" aria-controls="ai-engineering-zh" aria-pressed="false">中文</button>
  </div>

  <section id="ai-engineering-en" class="bilingual-post-panel is-active" data-bilingual-lang="en" lang="en" markdown="1">

Over the past month, I worked on a few curiosity-driven side projects that were not directly related to my computational biology or bioinformatics research. None of them were particularly complex, nor did they require an especially high technical barrier. But they made me realize something quite strongly: what I truly lack may not be mastery of a specific programming language, package, or command, but a more top-down engineering mindset.

I had seen the argument before: in the age of AI, the most important ability may no longer be writing code itself, but the ability to define and organize problems in an engineering-oriented way. For a long time, I did not feel this very deeply. To be honest, I also did not find AI dramatically more convenient for day-to-day bioinformatics analysis.

But when I was developing computational biology algorithms, I had already felt something similar. Code is only a small part of the work. The harder questions are: What problem am I actually trying to solve? How can this problem be decomposed into modules? Which parts should be built from scratch, and which parts can rely on existing tools? How should different modules communicate with each other? How do I know whether the entire system actually works, rather than merely looking successful in a local demo?

At the time, because AI was not yet powerful enough, these questions were partly obscured by the equally important challenge of writing the code itself. However, since late last year or early this year, the agent ecosystem has entered a phase of explosive growth. I suddenly realized that coding ability itself is becoming less central than before. This is especially true for computational biology, where our code usually does not need to deal with many real-world business constraints. In most cases, we mainly need to focus on the scientific question itself.

So, at least from the perspective of computational method development, engineering design has become a core issue that sits above the act of writing code. AI can help with implementation, but it is still very difficult for AI to replace me in defining the problem, structuring the system, and deciding what should be built.

Looking back, I think this may have something to do with my undergraduate training. A science background trains us very well to understand and explain a problem: why a phenomenon matters, what mechanism may underlie it, how we can test it through experiments or analysis, and how we eventually arrive at an explanatory answer.

But engineering problems follow a different logic. Engineering begins with a goal: Can I use the tools, code, models, and knowledge available to me to build a system that runs reliably? Scientific problems often pursue "why." Engineering problems first ask, "Can this be built, and can it work stably?" Only after that do we iterate.

A great example is [OpenClaw](https://github.com/openclaw/openclaw?utm_source=chatgpt.com). Today, OpenClaw and its related projects are still developing rapidly. But half a year ago, it was still just a relatively ordinary agent project based on [Pi](https://github.com/earendil-works/pi?utm_source=chatgpt.com), with the goal of connecting to messaging platforms. If we managed this kind of project purely with a scientist's mindset, perhaps it would never have been released at all. Every time it was about to ship, we might discover that its benchmark performance was not as strong as some newly released framework. Every time a new feature was added, we might also find it hard to fully explain why it worked. The project would then be postponed indefinitely.

Even when we are simply "building a tool," computational biologists may have very different design orientations. Some tools are driven by biological meaning, some by algorithmic novelty, and others by community extensibility. These differences may seem subtle, but even though I have not been in this field for that long, I can already feel them.

And these subtle differences in design orientation are precisely what begin to distinguish people when low-level code implementation becomes increasingly easy to delegate to AI. What truly determines the final product may no longer be how many syntactic details or command-line flags one can remember, but whether one can define a clear goal, design a reasonable structure, identify weak points in the system, and turn an idea into a genuinely usable tool through testing and iteration.

Starting from this post, I will share some reflections on engineering practices. The first topic in this series is SDD: **Spec-Driven Development**.

## Spec-Driven Development: From Slot-Machine Coding to Controlled System Building

One major reason many people resist using AI for serious work is that using AI often feels like gambling. You write a prompt, pull the lever, and never quite know what will come out. It may generate an impressive demo in one shot. But it may also choose the wrong tech stack, misunderstand your requirements, introduce strange abstractions, or turn what should have been a long-term maintainable project into a disposable script.

But this problem is not entirely caused by the model being "not smart enough." Very often, the real problem is that we have not given the model a stable, persistent, and verifiable project context.

GitHub's [Spec Kit](https://github.com/github/spec-kit?utm_source=chatgpt.com) describes SDD as an approach that puts specifications at the center of AI-assisted software development. The first sentence of their repository says:

> **An open source toolkit that allows you to focus on product scenarios and predictable outcomes instead of vibe coding every piece from scratch.**

Put simply, the idea explained in the GitHub Blog is very direct: a spec is no longer a document that is written once and then thrown away. Instead, it is a living artifact that evolves with the project. It is the shared source of truth among humans, tools, and AI agents. It records project details and development requirements, such as API conventions, approximate frontend UI and features, testing standards, and acceptance criteria.

Many people are familiar with Codex's Plan Mode. In my view, it can be understood as a one-off SPEC markdown file. A SPEC folder, by contrast, is more like a persistent collection of requirement documents and implementation guides (Table 1, [cite](https://www.augmentcode.com/guides/what-is-spec-driven-development)).

|Dimension|Vibe Coding|SDD|
|---|---|---|
|Primary artifact|Natural language prompts|Executable specifications|
|Scope|Full application generation|System-wide architectural contracts|
|Validation mechanism|Manual review, if any|Build fails on spec divergence|
|AI governance|None built in|Constitutional constraints and checkpoints|
|Where truth lives|Prompt history|Versioned specification|

## A Very Simple Example: Building a Web-Based Calculator

Let's use a very simple example. Suppose you are taking a class, and the assignment is to build a basic web-based calculator. Under the mindset of vibe coding, you might open an agent and directly type:

```
Please help me build a web-based calculator. I want it to be simple, but fully functional.
```

The agent will probably produce a calculator. But very soon, you may notice a series of problems: it uses JavaScript, while you prefer TypeScript; the color scheme is not what you had in mind; the buttons are arranged in a simple grid rather than laid out like a physical calculator; it does not consider mobile devices; and the state management design is not suitable for adding more operators later.

So you start adding requirements round after round, until the result finally looks good enough to submit.

The real problem appears next week, when the instructor asks you to continue building on the same calculator and add more operators and rules. You open the agent again, only to find that it has long forgotten the implicit agreements you went back and forth on last week. Maybe the color can still be reproduced, but the component structure, button layout, and extensibility of the calculation logic start drifting again.

Eventually, you realize that fixing the old project may cost more than rewriting it. For a small assignment, that may still be acceptable. But what if this were a large system?

This is the failure of prompt history as project memory.

SDD approaches the problem very differently. Even if you do not use Spec Kit, you can simply say to the agent:

```
I want to build a web-based calculator. Before writing any code, please ask me questions to define the project requirements. After I answer, write a detailed SPEC.md file as the long-lasting constitution for this project.
```

This step may look slow, but what it is really doing is turning vague requirements into stable context. A sufficiently capable agent will begin asking questions: Do you want TypeScript or JavaScript? Should React be used? Should the buttons imitate a physical calculator? Should keyboard input be supported? Do you need calculation history? Should parentheses and operator precedence be supported? How should the interface behave on mobile? How should invalid input be handled? Will scientific functions be added in the future?

You do not even need to know, at the very beginning, what aspects should be considered when building a calculator. You only need to answer the questions. The requirements gradually emerge through the conversation and are eventually organized into a SPEC file.

Next week, when the assignment asks you to continue developing the calculator, you ask the agent to update the SPEC first, and then generate the plan and tasks based on the updated SPEC. In this way, continuous development becomes possible.

## Why SDD Matters for Bioinformaticians

Bioinformatics is naturally compatible with SDD. The reason is simple: biological intent is already a form of specification.

The problem is that, in the past, these specifications were often scattered across many places. Some were in the biological question. Some were in the experimental design. Some were in the Methods section of a paper. Some were discussed during lab meetings. Some were written as comments in notebooks. And a large part of them existed only in our own heads.

For example, in a single-cell RNA-seq analysis project, the real requirement is never as simple as "help me run a Seurat pipeline." Behind that sentence, there are many questions:

```
- What biological question is this project actually trying to answer?
- Which conditions, time points, tissues, or treatments do the samples come from?
- Which samples should be included, and which should be excluded?
- How should QC thresholds be set? Why?
- Should doublet detection be performed?
- How should batch effects be handled?
- How should the clustering resolution be chosen?
- What criteria should define marker genes?
- What are the comparison groups for differential expression?
- Which intermediate results must be saved?
- Should the final deliverable be figures, tables, a report, or a reusable pipeline?
- What counts as success? What indicates that the pipeline has failed?
```

For bioinformatics projects, the most dangerous errors are often not cases where the code crashes. The most dangerous errors are cases where the code runs without any error, while the analysis goal has quietly drifted.

When people talk about an "experienced bioinformatics engineer," I think this is exactly what they mean: someone who can complete the project robustly without letting the goal drift. By contrast, the biological interpretation of the results may actually be something that scientists are better positioned to handle.

Therefore, I think the significance of SDD for bioinformaticians is not merely that it makes AI-generated code more stable. More importantly, it helps us translate biological intent into engineering constraints.

It allows us to persistently record the biological meaning of a project, the overall pipeline design, parameter choices, QC standards, replaceable modules, intermediate results, and acceptance criteria. In this way, the agent is no longer just a temporary assistant that helps you write scripts. It becomes more like a collaborator working under the same project constitution.

Of course, learning how to write good SPEC files also requires experience and foundational knowledge. In fact, the entry barrier may be no lower than when we first learned to write code by hand. But I believe this working logic is absolutely part of the future.

For computational biology algorithm development, some tasks are even more naturally suited to this approach. Benchmarking is one example: the goal is clear, the tools are clear, the metrics are clear, and the acceptance criteria are clear. In principle, after receiving a good SPEC, an agent can complete all the experiments. I would even say that, in some cases, its stability and objectivity may be better than a human's.

In short, we should open up our thinking. For tasks that are suitable for agents, we should learn to assign them through a SPEC system.

As a famous Chinese saying goes: taking time to sharpen the axe does not delay the work of chopping wood. In the age of AI coding, perhaps this saying can offer us a new kind of inspiration.

  </section>

  <section id="ai-engineering-zh" class="bilingual-post-panel" data-bilingual-lang="zh" lang="zh-CN" hidden markdown="1">

这一个月里，我做了一些和计算生物、生物信息科研任务并不直接相关的兴趣项目。它们不一定复杂，也不一定有多高的技术门槛，但它们让我很强烈地意识到：我真正欠缺的，可能并不是某一种编程语言、某一个包、某一条命令，而是一种更自上而下的工程思维。

以前我对"AI 时代最重要的能力不再是写代码本身，而是工程化地定义和组织问题"这个观点并没有太深的体会，其实我也没有发现在解决生物信息分析的方面有多方便。回过头看，我在开发计算生物学算法的时候，很早就有所体会，代码只是很小的一部分，更困难的是：我到底要解决什么问题？这个问题可以被拆成哪些模块？哪些部分可以调用现成工具？不同模块之间如何衔接？如何判断整个系统是真的 work，而不是只是在一个局部 demo 里看起来 work？这些问题，当时由于AI并没有很强大，所以被编码能力的同等重要性给掩盖了。然而，从今年或者去年年末开始，Agent迎来了井喷的爆发阶段，我突然发现代码能力本身不再是唯一的核心壁垒，尤其是我们计算生物学的代码，并不需要考虑各种现实业务因素，只要足够关注科学问题即可。所以单纯从计算方法编程角度，工程学设计问题，成为了凌驾于代码编写能力之上的一个核心问题：AI 可以辅助，但很难替我完成。

回头看，我意识到这也许和我的本科训练有关。Science 背景的训练更擅长教我们如何理解和解释一个问题：为什么这个现象重要，它背后的机制可能是什么，我们如何通过实验或分析验证它，最后如何给出一个解释性的答案。但它较少训练我们如何从一个目标出发，把已有工具、代码、模型和知识组织成一个稳定运行的系统。科学问题往往追求"为什么"，而工程问题首先追求"能不能构建出来，并且稳定地工作"，然后再慢慢迭代。一个极佳的例子就是[OpenClaw](https://github.com/openclaw/openclaw)，直到今天它和它的衍生项目仍然在高速发展，而它半年前，还仅仅是基于[pi](https://github.com/earendil-works/pi)的一个目的是接入聊天软件的普通agent项目：如果我们用科学家的思维方式管理，那么也许这个项目永远不会出现：每次将要发出来的时候，都会发现自己benchmark跑不过新发的框架，然后每次新加的功能为什么work也无法很好解释；所以就被无限期的推迟了。

同样是开发一个工具，计算生物学家开发的工具之中，有生物学意义导向的，有算法学导向的，也有社区拓展性导向的......这些差异看似微妙，但我即便混迹于此没有特别久，也多少有体会。这些微妙的设计导向区别，就是当底层代码实现变得越来越容易被 AI 辅助完成时，真正区分人的能力。真正决定最终成品的也许就不再是记住多少语法和命令，而是能否提出清晰的目标、设计合理的结构、识别系统的薄弱环节，并在不断测试和迭代中把一个想法变成真正可用的工具。从这一期博客开始，我也将分享一些工程学实践的感悟。这个系列的第一个分享点，就是 SDD：Spec-Driven Development。

## Spec-Driven Development：从"抽卡式写代码"到"可控地构建系统"

很多人排斥 AI 用于正式工作的一个核心原因是：每次使用 AI 都像在赌博。你写下一段 prompt，拉动拉杆，然后永远不知道吐出来的会是什么。它可能一次性写出一个惊艳的 demo，也可能在你没注意的地方选错技术栈、误解需求、引入奇怪的抽象，或者把一个本来应该长期维护的项目写成一次性脚本。但这个问题不完全是模型不够聪明。很多时候，是我们没有给它一个稳定、持久、可验证的项目上下文。

GitHub的[SpecKit](https://github.com/github/spec-kit)把SDD描述为一种把 specification 放在 AI-assisted software development 中心的方法。 这是来自他们仓库的第一句话：

> **An open source toolkit that allows you to focus on product scenarios and predictable outcomes instead of vibe coding every piece from scratch.**

简而言之，GitHub Blog 对这个思想的解释也很直接：spec 不再是写完就扔掉的文档，而是一个会随着项目演化的 living artifact，是人、工具和 AI agent 共享的 source of truth，是一种记录了项目细节和开发要求的文档，比如约定API，约定前端的大致UI和feature，测试验收标准，等等。很多人知道codex 的plan mode，其实本质上就是单次的SPEC markdown；而SPEC folder，更像是持久的一系列需求文档和实现指南（Table 1， [cite](https://www.augmentcode.com/guides/what-is-spec-driven-development))。

| Dimension            | Vibe Coding                 | SDD                                        |
| -------------------- | --------------------------- | ------------------------------------------ |
| Primary artifact     | Natural language prompts    | Executable specifications                  |
| Scope                | Full application generation | System-wide architectural contracts        |
| Validation mechanism | Manual review (if any)      | Build fails on spec divergence             |
| AI governance        | None built-in               | Constitutional constraints and checkpoints |
| Where truth lives    | Prompt history              | Versioned specification                    |

## 非常简单的例子：开发一个网页计算器

这里举一个很简单的例子，你选了一门课，作业是一个最简单的网页计算器；在vibe coding的认知下，也许就会直接打开一个agent：

```markdown
> Please help me to build a web-based calculator. I want it to be simple but with full functions.
```

Agent 很可能真的会做出一个计算器。但你很快发现：它用了 JavaScript，而你更想用 TypeScript；它的颜色不是你想要的；按钮排布只是普通 grid，而不是模拟实体计算器；它没有考虑移动端；状态管理也不利于后续加入更多运算规则。于是你开始一轮又一轮地补充需求，直到它终于看起来能交作业。

问题是，下周老师要求你在这个计算器基础上继续开发，加入更多运算符和规则。你再次打开 agent，却发现它早就不记得你们上周那些来回修改的隐含约定。颜色也许还能复现，但组件结构、按钮布局、运算逻辑扩展方式又开始漂移。最后你发现，与其继续修，不如重写一个：所以导致，维护成本反而超过了重写成本。这对于小项目来说还算可以接受，但是一个大系统项目呢？这就是 prompt history 作为项目记忆的失败。

但是，如果你知道了SPEC driven development，即便不用spec-kit，而 SDD 的做法会完全不同，只需要 agent 说：

```
I want to build a web-based calculator. Before writing any code, please ask me questions to define the project requirements. After I answer, write a detailed SPEC.md file as the long-lasting constitution for this project.
```

这一步看似很慢，但它实际上是在把模糊需求变成稳定上下文。一个足够好的 agent 会开始问你问题：你要 TypeScript 还是 JavaScript？是否使用 React？按钮是否需要模拟实体计算器？是否支持键盘输入？是否需要历史记录？是否要支持括号和运算优先级？移动端怎么显示？错误输入如何处理？未来是否会加入科学计算函数？

你甚至不需要一开始就知道"开发一个计算器应该考虑哪些方面"。你只需要回答问题。需求会在问答中逐渐浮现，最后被整理成一个 SPEC 文件。到了下周，当作业要求你继续开发这个计算器时，你就让 agent 先更新 SPEC，再基于更新后的 SPEC 生成 plan 和 tasks。这样，持续开发才成为可能。

## 为什么SDD对生物信息学家是重要的

生物信息学天然和SDD是契合的。原因很简单：生物学意义需求本身就是一套 specification。只是过去这些 specification 往往散落在不同地方：一部分在 biological question 里，一部分在实验设计里，一部分在 paper 的 methods 里，一部分在 lab meeting 的讨论里，一部分在 notebook 的注释里，还有很大一部分只存在于我们自己的脑子里。

比如一个单细胞 RNA-seq 分析项目，真正的需求从来不是"帮我跑一个 Seurat pipeline"这么简单。它背后至少包含：

```
- 这个项目真正想回答的 biological question 是什么？
- 样本来自哪些条件、时间点、组织或处理组？
- 哪些样本应该被纳入，哪些应该被排除？
- QC 阈值如何设定？为什么？
- doublet detection 是否需要做？
- batch effect 应该如何处理？
- clustering resolution 如何选择？
- marker gene 的判定标准是什么？
- differential expression 的比较组是什么？
- 哪些中间结果必须保存？
- 最终交付物是 figures、tables、report，还是可复用 pipeline？
- 什么结果算成功？什么结果说明 pipeline 失败？
```

对生物信息项目来说，最危险的错误往往不是代码报错，而是代码没有报错，但分析目标已经悄悄漂移了。"有经验的生信工程师"说的就是能在目标不发生漂移的情况下，稳健完成这个项目任务；生物学的需求和结果解读其实反而可能是科学家更能胜任的。所以我认为，SDD 对生物信息学家的意义不只是"让 AI 写代码更稳定"，而是帮助我们把 biological intent 转化成 engineering constraints。它可以把项目的生物学意义、整体 pipeline 设计、参数选择、QC 标准、可替换模块、中间结果和验收标准持久化下来。这样 agent 不再只是一个临时帮你写脚本的助手，而更像是一个围绕同一套项目宪法工作的合作者。当然，这个学习如何编写SPEC的过程，也是需要比较丰富的经验和基础的，甚至入门难度不亚于最初我们手写代码，但是我认为这套工作逻辑绝对是未来的发展趋势。

那么计算生物学算法开发，有些任务就更适合了，比如benchmark：目的明确，工具明确，metrics明确，验收标准明确，所有的实验理论上都可以由agent接收SPEC后完成：我甚至觉得稳定性和客观性比人还好...

总之，应该要打开思路，把适合agent完成的任务，用SPEC system的形式下发，中国古话言"磨刀不误砍柴工"，在AI coding时代也许能给各位新的启发。

  </section>
</div>

<style>
.bilingual-post-switch {
  display: flex;
  flex-wrap: wrap;
  gap: 0.5rem;
  margin: 0 0 1.25rem;
  padding: 0.5rem;
  border: 1px solid var(--global-divider-color);
  border-radius: 8px;
  background: var(--global-bg-color);
}
.bilingual-post-button {
  border: 1px solid var(--global-divider-color);
  border-radius: 999px;
  padding: 0.35rem 0.8rem;
  background: transparent;
  color: var(--global-text-color);
  font-size: 0.86rem;
  font-weight: 700;
  cursor: pointer;
}
.bilingual-post-button.is-active {
  border-color: var(--global-theme-color);
  background: var(--global-theme-color);
  color: #fff;
}
.bilingual-post-panel[hidden] {
  display: none !important;
}
.bilingual-post-panel {
  min-width: 0;
}
.bilingual-post-panel img {
  display: block;
  width: 100%;
  max-width: min(100%, 760px);
  height: min(52vh, 460px);
  object-fit: contain;
  object-position: center;
  margin: 1.25rem auto;
  border-radius: 8px;
}
@media (max-width: 640px) {
  .bilingual-post-panel img {
    max-width: 100%;
    height: min(42vh, 360px);
  }
}
</style>

<script>
(() => {
  document.querySelectorAll('.bilingual-post').forEach((root) => {
    const buttons = root.querySelectorAll('[data-bilingual-target]');
    const panels = root.querySelectorAll('[data-bilingual-lang]');
    const title = root.closest('.post')?.querySelector('.post-title');

    buttons.forEach((button) => {
      button.addEventListener('click', () => {
        const target = button.dataset.bilingualTarget;

        buttons.forEach((candidate) => {
          const active = candidate === button;
          candidate.classList.toggle('is-active', active);
          candidate.setAttribute('aria-pressed', active ? 'true' : 'false');
        });

        panels.forEach((panel) => {
          const active = panel.dataset.bilingualLang === target;
          panel.hidden = !active;
          panel.classList.toggle('is-active', active);
        });

        if (title) {
          title.textContent = target === 'zh' ? root.dataset.titleZh : root.dataset.titleEn;
        }
      });
    });
  });
})();
</script>
