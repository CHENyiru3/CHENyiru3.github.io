---
layout: post
title: "How the Research Radar Works: An Automated Literature Curation Pipeline (Not Just for Fun)"
date: 2026-05-15
description: "A behind-the-scenes look at the automated pipeline that scans academic RSS feeds, curates recommendations with AI, and renders a daily research digest — from cron jobs to Jekyll templates."
tags:
  - tech
  - automation
  - research-radar
categories:
  - tech-blog
bilingual: true
default_lang: en
related_posts: false
title_zh: "Research Radar 工作原理：一套自动化学术文献策展管线（不止是玩玩而已）"
description_zh: "幕后揭秘：一套自动化管线如何扫描学术 RSS 源、用 AI 策展推荐论文、并渲染出每日研究简报——从 cron 任务到 Jekyll 模板。"
---

<div
  class="bilingual-post"
  data-default-lang="en"
  data-title-en="How the Research Radar Works: An Automated Literature Curation Pipeline (Not Just for Fun)"
  data-title-zh="Research Radar 工作原理：一套自动化学术文献策展管线（不止是玩玩而已）">
  <div class="bilingual-post-switch" role="group" aria-label="Language switcher">
    <button class="bilingual-post-button is-active" type="button" data-bilingual-target="en" aria-controls="research-radar-en" aria-pressed="true">English</button>
    <button class="bilingual-post-button" type="button" data-bilingual-target="zh" aria-controls="research-radar-zh" aria-pressed="false">中文</button>
  </div>

  <section id="research-radar-en" class="bilingual-post-panel is-active" data-bilingual-lang="en" lang="en" markdown="1">

Every morning at 09:15, a scheduled task starts running on my computer.

It scans the local NetNewsWire database for new papers from the past 24 hours, pulls unread entries from more than twenty journals and preprint sources, then sends the titles and abstracts to a model so it can help me choose the 13 to 15 papers most worth reading that day. I call this small agent Clawdie. The name comes from Claude Code, because that was the first tool that really pulled me into the AI agent rabbit hole.

The scheduling layer runs on Hermes Agent, the timing is controlled by cron, and the summarization and judgment are done by DeepSeek-V4-Pro. The final result is written into my Jekyll website repository, then GitHub Pages builds it automatically. By the time I open my computer in the morning, the Research Radar page is usually already updated.

It does not read papers for me. Its job is narrower: remove noise, filter out papers that are completely unrelated to my interests, and leave me with the day's themes, recommendation reasons, and DOI links. My real job is still to decide which papers deserve opening and reading carefully.

## Why I Built This

Academic publishing really does move too fast.

Nature, Science, Cell, and all of their sibling journals publish new articles every day. Add bioRxiv, medRxiv, and specialized journals on top of that, and manually checking journal websites quickly becomes mechanical labor. The return is not even that good.

The hardest part is not simply that there are too many papers. It is that most papers have nothing to do with the question I actually care about on a given day. I tried several approaches before: RSS readers, email alerts, Zotero keyword searches, and posts people shared on social media. I still use RSS, mainly NetNewsWire. It works well, but the problem is straightforward: I still have to screen everything myself.

![The RSS reader I use every day on macOS: NetNewsWire.](https://i.imgur.com/lAzG0vx.png)

Looking at dozens of titles every day, opening some of them, scanning abstracts, and deciding "irrelevant", "maybe relevant", or "later" is not hard in isolation. Repeating it every day consumes attention. Worse, the few truly important papers are often mixed into a pile of unrelated material, which makes them easy to miss.

So the goal of Research Radar is simple. It automates the most annoying first layer of screening: the machine reads the abstracts first and ranks them by my research interests; I decide which ones are worth serious reading. I do not think an agent can tell you which scientific questions are important, so it can only serve as an initial filter, helping you spend time on more important problems.

## What the Daily Page Looks Like

The morning digest is roughly split into two parts.

The first part is the day's hot topic. Clawdie looks at the papers selected that day, combines them with my recent Zotero additions and the research-direction hints written into the system prompt, and finds the clearest thread connecting several papers. It usually generates a short summary and a few signal points.

This part is useful to me because it is not just a paper list. It answers a more practical question: among today's new papers, is there any trend worth noticing?

The second part is the article list. Right now I split it into several sections:

- **Computational**, usually 5 papers, mainly methods, algorithms, models, and AI tools.
- **Biomedicine**, usually 5 papers, focused on biological findings related to my research.
- **Other Fields**, usually 3 papers, mostly AI or computational work that does not fully belong to biomedicine but still seems worth a glance.

On Fridays, there is an extra Biotech News Delivery section with about 5 items, mainly industry updates.

Each paper includes several fields: a very short contribution summary, a "why it matters" note, a more personal "why for Yiru" note, and a recommendation level, such as read carefully, skim, or awareness.

![Research Radar example: a daily digest with hot topic and curated article sections.](https://i.imgur.com/FfNakde.png)

I deliberately ask the model to write "why", rather than just give a score. A plain relevance score is not that helpful, because it is hard for me to make a real decision from the difference between 8 and 9. But if it explains how a paper connects to spatial omics, single-cell analysis, computational immunology, or the tumor microenvironment, I can decide much faster whether to open the original paper.

## The Workflow

The system is not especially complicated. At a high level, it has three parts: collect sources, let the model curate, and generate the web page.

### 1. Sources: Control the Entry Point First

All paper sources enter NetNewsWire first.

I subscribe to more than twenty RSS feeds, mainly including AOP feeds from the Nature family, Cell Press, Science, and several computational biology journals and preprint sources, such as PLOS Computational Biology and the methods-oriented sections of bioRxiv.

These feeds are organized into folders in NetNewsWire. The advantage is that I do not need to maintain a separate paper crawling system. NetNewsWire already handles RSS updates, read/unread state, and local storage for me. More conveniently, the NetNewsWire article database is just a local SQLite file. This is important. My pipeline can read that SQLite database directly, without applying for API keys, using OAuth, or worrying that some third-party service will change its interface one day. As long as NetNewsWire updates normally, Research Radar can get the latest article state.

In short, NetNewsWire is responsible for "what is new"; my script is responsible for "what is worth reading".

### 2. AI in the Loop: Let the Model Read Abstracts First

Every day at 09:15, a Hermes Agent cron task starts. It connects to the NetNewsWire SQLite database, retrieves unread articles from the past 24 hours, formats their titles, sources, abstracts, and links as input, and sends them to DeepSeek-V4-Pro.

The prompt contains my research profile. Right now it mainly includes spatial omics, single-cell analysis, computational immunology, tumor microenvironment, biomedical AI, and the use of large models and foundation models in scientific research.

The model does more than score papers. It needs to complete several steps:

- First, judge how relevant each paper is to my research interests.
- Then, select the papers most worth reading for each section.
- Next, write a contribution summary, importance explanation, and personal relevance note for each paper.
- Finally, extract the day's hot topic from the selected papers.

These contents are written into a Markdown file with YAML metadata at the top. This format works well for later Jekyll rendering.

There is one practical problem, though: model-generated YAML is not always stable.

Early on, it once missed a closing `---`, and the entire digest rendered as a blank page on the website. Sometimes the array format was wrong, or the DOI field was empty. So I later added a second cron task that runs at 10:00 specifically for QA.

The QA task checks YAML structure, delimiters, required fields, DOI values, array format, and similar details. If the issue is small, it fixes it automatically. If the issue is more serious, the bad file should not reach the live page.

It is usually invisible, but necessary. After a week, it has not caused any problems.

### 3. Publishing: From Markdown to the Website

The Markdown file generated by the model is placed in the website repository under `_research_radar/YYYY-MM-DD.md`.

Hermes Agent commits this file and pushes it to GitHub. After that, GitHub Pages automatically rebuilds the website. Under normal conditions, it takes only a few dozen seconds from push to live page.

On the Jekyll side, I use collections to manage Research Radar entries, then Liquid templates generate the index page, individual pages, and RSS feed.

The final result is a fixed-format digest every day, automatically archived, with permanent links and RSS subscription support.

I try to keep this part boring. Content generation is already unstable enough; the publishing path does not need extra tricks.

## Design Tradeoffs

This system has been able to run stably mostly because of several conservative choices.

The first is using SQLite instead of connecting a pile of APIs.

NetNewsWire's local database already records which articles are new and which ones have been read. Reading it directly is enough. There are no API rate limits, no token rotation, and no sudden breakage because some service redesigned its interface. For this kind of personal automation, a local file is actually the most reliable interface.

The second is asking the model to summarize, not just filter.

At first I thought it would be enough to ask the model to decide "relevant" or "irrelevant". Later I realized that was not enough. What I really needed was not a binary classifier, but a filter that could explain its judgment.

I want to know why a paper is worth reading, or why it can be skipped. If the explanation is clear enough, then even when I decide not to read the paper, I still roughly know what I am missing.

The third is separating generation, checking, and display.

The model only generates content. The QA task only checks structure and repairs problems. Jekyll only renders the content into pages. That means any layer can be replaced independently.

For example, if I later want to switch from DeepSeek to another model, the website template does not need to change. If I want to redesign the Research Radar page, the prompt does not need to change.

The fourth is not letting automation cross the boundary.

Research Radar only adds papers, ranks them, labels them, and writes explanations. It does not automatically add papers to Zotero, decide research directions for me, or delete anything.

I still want the final judgment to stay with me. Automation can save attention, but it should not use attention on my behalf.

## How It Feels So Far

Research Radar is not a complicated system. In essence, it is cron, SQLite, one model call, Markdown, Jekyll, and GitHub Pages stitched together.

But it solves a problem I face every day: I do not want to spend the clearest part of my morning filtering out dozens of irrelevant papers.

Now that part goes to the machine first. It narrows the scope, gives reasons, and puts the result somewhere I already check every day. Which paper to actually read, how to read it, and whether to follow up afterward are still my decisions.

I think this is the right place for this kind of research automation. It is not about letting the machine do research for you, but about letting it waste a little less of your attention. I hope this post helps, and I am happy to talk more about it. This pipeline is not open source for now, because it feels fairly simple, and the idea here should be enough to build something similar.

  </section>

  <section id="research-radar-zh" class="bilingual-post-panel" data-bilingual-lang="zh" lang="zh-CN" hidden markdown="1">

每天早上 09:15，我的电脑上会有一个定时任务开始跑。

它会去 NetNewsWire 的本地数据库里扫一遍过去 24 小时的新论文，从二十多个期刊和预印本源里抓出未读条目，然后把标题和摘要交给模型，让它帮我挑出当天最值得看的 13 到 15 篇。这个小代理我叫它 Clawdie，名字来自 Claude Code，因为它算是第一个把我带进 AI agent 坑里的工具。

调度这部分跑在 Hermes Agent 上，时间由 cron 控制，真正做摘要和判断的是 DeepSeek-V4-Pro。最后生成的结果会被写进我的 Jekyll 网站仓库，GitHub Pages 自动构建。等我早上打开电脑时，Research Radar 页面通常已经更新好了。

它不会替我读论文。它做的事情更窄一点：去掉噪声，过滤掉和我的关注点完全不相关的，我真正需要做的，是扫一眼当天的主题、推荐理由和 DOI 链接，然后决定哪些论文值得打开原文细读。

## 为什么要做这个

学术出版的速度确实太快了。

Nature、Science、Cell 以及它们的一堆子刊，每天都会冒出不少新文章。再加上 bioRxiv、medRxiv 和各种专业期刊，光靠人手刷这些期刊官网，很快就会变成一件机械劳动，而且成效很难说很好。

最麻烦的地方不是“论文太多”。而是绝大多数论文和我当天真正关心的问题没有关系。我以前试过几种办法：RSS 阅读器、邮件提醒、Zotero 关键词搜索，还有社交媒体上的转发。RSS 我现在也还在用，主要是 NetNewsWire。它很好用，但问题也很直接：最后还是要我自己一篇篇筛。

![我每天在 macOS 上用的 RSS 阅读器——NetNewsWire。](https://i.imgur.com/lAzG0vx.png)

每天面对几十篇标题，点开、扫摘要、判断“不相关”“可能相关”“以后再说”，这件事单独看不难，但每天重复就很消耗注意力。更糟的是，真正重要的几篇文章经常混在一堆无关内容里，很容易被漏掉。

所以 Research Radar 的目标很简单，它只是把每天最烦人的第一层筛选自动化掉：机器先读一遍摘要，按我的研究兴趣排个序；我再决定哪些值得认真读。我不认为Agent能帮你思考什么科研问题是重要的，所以它只能充当一个最初的筛选器，帮助你把时间放在更重要的问题上

## 每天生成的页面长什么样

每天早上的简报大概分成两部分。

第一部分是当天的热点主题。Clawdie 会看当天被选中的论文，结合我最近的 Zotero 新增条目和系统里写好的研究方向提示，找出几篇论文之间最明显的主线。它通常会生成一小段总结，再列几个信号点。

这个部分对我来说很有用，因为它不是简单列论文，而是回答一个更现实的问题：今天这些新文章里，有没有什么趋势值得注意？

第二部分是文章列表。现在我把它分成几个版块：

- **计算方向**，通常 5 篇，主要是方法、算法、模型和 AI 工具。
- **生物医学方向**，通常 5 篇，偏向和我研究相关的生物学发现。
- **其他领域**，通常 3 篇，主要放一些不完全属于生物医学、但我觉得可能值得瞄一眼的 AI 或计算方向文章。

周五会多一个 Biotech 新闻速递，大概 5 条，主要是行业动态。

每篇文章下面会有几类信息：一段很短的贡献总结，一段“为什么重要”，一段更个人化的“为什么对 Yiru 重要”，还有一个推荐等级，比如精读、浏览、知晓。

![Research Radar 示例——一份包含热点话题和策展文章分区的每日简报。](https://i.imgur.com/FfNakde.png)

这里我刻意让模型写“为什么”，而不是只给分数。单纯的相关性分数其实帮助不大，因为所谓 8 分和 9 分相关性之间我很难形成判断。但如果它能说明这篇文章和空间组学、单细胞分析、计算免疫学或者肿瘤微环境之间的关系，我就能更快决定要不要点进去读原文。

## 整个流程

这个系统没有特别复杂。大体上就是三段：收集来源、模型策展、生成网页。

### 1. 信息源：先把入口管好

所有论文来源先进入 NetNewsWire。

我订阅了二十多个 RSS 源，主要包括 Nature 系列期刊的 AOP、Cell Press、Science，以及一些计算生物学相关的期刊和预印本源，比如 PLOS Computational Biology、bioRxiv 的方法方向等。

这些源在 NetNewsWire 里按文件夹组织。这样做的好处是，我不需要另外维护一套文章抓取系统。NetNewsWire 本身已经帮我处理了 RSS 更新、已读未读状态和本地存储。更方便的是，NetNewsWire 的文章数据库就是一个本地 SQLite 文件。这点很关键。因为我的管线可以直接读这个 SQLite 数据库，不需要申请 API key，不需要 OAuth，也不用担心某个第三方服务哪天改接口。只要 NetNewsWire 正常更新，Research Radar 就能拿到最新的文章状态。

简单说，NetNewsWire 负责“有什么新文章”，我的脚本负责“哪些值得看”。

### 2. AI介入：让模型先读一轮摘要

每天 09:15，Hermes Agent 的 cron 任务会启动。它会连接 NetNewsWire 的 SQLite 数据库，取出过去 24 小时的未读文章，然后把标题、来源、摘要和链接整理成输入，交给 DeepSeek-V4-Pro。

提示词里写了我的研究画像。现在主要包括这些方向：空间组学、单细胞分析、计算免疫学、肿瘤微环境、生物医学 AI，以及大模型和基础模型在科研里的应用。

模型要做的事情不只是打分。它需要完成几步：

- 先判断每篇文章和我的研究兴趣有多相关。
- 再按不同版块挑出当天最值得看的文章。
- 然后为每篇文章写贡献摘要、重要性解释和个人相关性说明。
- 最后再从被选中的文章里提炼当天的热点主题。

这些内容会被写成一个 Markdown 文件，文件头部带 YAML metadata。这个格式很适合 Jekyll 后续渲染。

不过这里也有一个很实际的问题：模型生成 YAML 并不总是稳定。

早期有一次，它少写了一个闭合的 `---`，结果整份简报在网站上直接变成空白。还有时候数组格式不对，或者 DOI 字段为空。于是我后来加了第二个 cron 任务，10:00 跑一遍，专门做 QA。

这个 QA 任务会检查 YAML 结构、分隔符、必填字段、 DOI、数组格式等。如果只是小问题，它会自动修掉；如果问题比较严重，就不要让坏文件进入线上页面。

它平时没有存在感，但很有必要。如今过了一周了，没有出现任何问题

### 3. 发布：从 Markdown 到网页

模型生成的 Markdown 文件会被放到网站仓库的 `_research_radar/YYYY-MM-DD.md` 里。

Hermes Agent 会把这个文件 commit，然后 push 到 GitHub。之后 GitHub Pages 会自动重新构建网站。正常情况下，从 push 到页面上线大概几十秒。

Jekyll 这边我用了 collections 来管理 Research Radar 的条目，再用 Liquid 模板生成索引页、单篇页面和 RSS feed。

所以最后的效果是：每天一份固定格式的简报，自动归档，有永久链接，也可以继续被 RSS 订阅。

这部分我尽量做得无聊一点。内容生成已经足够不稳定了，发布链路就不要再花里胡哨。

## 几个设计取舍

这个系统能稳定跑起来，主要靠几个比较保守的选择。

第一个是用 SQLite，而不是再接一堆 API。

NetNewsWire 的本地数据库已经记录了哪些文章是新的、哪些已经读过。直接读它就够了。没有 API 限流，没有 token 轮换，也不会因为某个服务改版而突然坏掉。对这种个人自动化来说，本地文件反而是最可靠的接口。

第二个是让模型做总结，而不是只做过滤。

一开始我以为只要让模型判断“相关 / 不相关”就行。后来发现这不够。因为我真正需要的不是一个二分类器，而是一个能解释判断过程的筛选器。

我想知道一篇论文为什么值得看，或者为什么可以跳过。只要解释写得足够清楚，就算我最后不读那篇文章，我也大概知道自己错过的是什么。

第三个是把生成、检查和展示分开。

模型只负责生成内容。QA 任务只负责检查结构和补救。Jekyll 只负责把内容渲染成页面。这样任何一层都可以单独替换。

比如我以后想把 DeepSeek 换成别的模型，网站模板不用动。或者我想重新设计 Research Radar 页面，提示词也不用改。

第四个是不要让自动化越界。

Research Radar 只会添加文章、排序、打标签和写解释。它不会自动把论文加入 Zotero，不会替我决定研究方向，也不会删除什么东西。

我还是想把最后的判断留给自己。自动化可以帮我省注意力，但不应该替我使用注意力。

## 目前的感受

Research Radar 不是一个复杂系统。它本质上就是 cron、SQLite、一次模型调用、Markdown、Jekyll 和 GitHub Pages 拼起来的东西。

但它解决了一个我每天都会遇到的问题：我不想把早上最清醒的时间花在筛掉几十篇无关论文上。

现在这部分工作先交给机器。它帮我把范围缩小，给出理由，再把结果放到一个我每天都会看的地方。至于真正读哪篇、怎么读、读完之后要不要追下去，还是我自己决定。

我觉得这才是这类科研自动化比较合适的位置。不是让机器替你做研究，而是让它少浪费一点你的注意力。希望这篇博客能帮到你，也欢迎任何交流！目前这个pipeline没有开源，因为感觉似乎比较简单，按这个思路都可以完成。

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
