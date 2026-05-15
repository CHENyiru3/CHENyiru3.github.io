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

Hey there, blog readers! 👋 Welcome to this corner of the internet. Glad you stopped by.

Every morning at 09:15, while I'm still on my first coffee, an AI agent I named **Clawdie** (named after Claude Code, the first agent that got me into this stuff) scans freshly published papers across two dozen journals, reads their abstracts, and decides which 13 to 15 are worth my attention. Clawdie runs on [Hermes Agent](https://github.com/nousresearch/hermes-agent) as the brain that schedules everything, with cron jobs handling the timing, and it hands the heavy reading and reasoning to **DeepSeek-V4-Pro** (love DS for this kind of work).

By the time I open my laptop, a ranked digest is waiting on my [Research Radar page](https://chenyiru3.github.io/research-radar/) — summaries, relevance explanations, and a "hot topic" that pulls together the day's biggest theme. I scan the lists and recommendation reasons, and click through to the actual DOI link when an article looks like it matters to my research or interests. Attention is all I need, and Research Radar saves me some time.

## Why I Built This

Academic publishing moves fast. Across Nature, Science, Cell, and their sibling journals alone, papers appear every day — and that's before you count preprints. You can't read them all. Most are irrelevant to my work.

But sometimes, a handful are essential. When the critical information is buried by the noise, missing it would genuinely suck.

![The RSS reader I use every day on macOS — NetNewsWire.](https://i.imgur.com/lAzG0vx.png)

I tried RSS readers (I use [NetNewsWire](https://github.com/Ranchero-Software/NetNewsWire/), by the way), email alerts, and Zotero keyword searches. Even social media — Twitter, WeChat, Xiaohongshu — but you know how that goes. They all had the same problem: they required *me* to triage. Every. Single. Day. Going through 50+ papers deciding "not relevant, not relevant, not relevant" just to find the 3 that matter is exhausting.

The Research Radar flips this — the machine does the triage, I do the reading.

## What the Daily Digest Looks Like

Before I get into how it works, here's what lands on the page each morning.

**Hot Topic** — Clawdie picks up on my recent reading patterns and Zotero additions, combined with a system-level prompt about my research areas, and pulls together the day's strongest thread connecting 2–5 of the selected papers. It gives me a short narrative plus a few signal bullet points — basically a one-paragraph answer to "what should I be paying attention to today?"

**Curated Article Sections** — The digest breaks down into three categories (plus a Friday bonus):

- **Computational** (5 articles) — methods, algorithms, and AI tools
- **Biomedicine** (5 articles) — biological discoveries relevant to my research
- **Other Fields** (3 articles) — AI breakthroughs from outside biomedicine worth knowing about
- **BioTech News Delivery** (5 articles, Fridays only) — industry news for context

For each article, I get a concise summary, a "Why it matters" field-level explanation, a personalized "Why for Yiru" note, and a recommendation level (READ FULL / SKIM / AWARENESS). If something catches my eye, I click the DOI link and read the real paper.

![Research Radar example — a daily digest with hot topic and curated article sections.](https://i.imgur.com/FfNakde.png)

## The Pipeline: How to Reproduce It

Three steps. Here's how you can set up your own.

### Step 1: Feed Collection — Control Your Source Quality

Curated RSS feeds from 20+ journals and preprint servers flow into [NetNewsWire](https://netnews.wire.com/), a native macOS RSS reader. I organize the feeds into folders: AOP (Accepted Online Papers) for Nature journals, Cell Press journals, Science, plus computational biology venues like PLOS Computational Biology and bioRxiv's methods section.

NetNewsWire stores its article database as a local SQLite file. Why does this matter? Because the entire feed state is queryable without any API keys, OAuth tokens, or third-party services. No accounts to manage, nothing to rotate. Just a SQLite database on disk.

### Step 2: AI Curation — Let the Model Read and Sort

At 09:15, a [Hermes Agent](https://hermes-agent.nousresearch.com/) cron job wakes up, connects to the NetNewsWire SQLite database, grabs all unread articles from the last 24 hours, and feeds them to DeepSeek-V4-Pro with a detailed curation prompt.

The prompt encodes my research profile — right now: spatial omics, single-cell analysis, computational immunology, tumor microenvironment biology, biomedical AI, and foundation models. The model scores each article for relevance, picks the top ones for each section, and generates all the summaries, explanations, hot topic synthesis, and recommendation labels you saw above.

Behind the scenes, a second cron job fires at 10:00 as a QA watchdog. AI-generated YAML is fragile — early on, a missing closing `---` once turned an entire digest into a blank page. The watchdog checks the YAML structure (no missing delimiters, no null DOIs, proper array formatting) and auto-patches any issues before they hit the live site. When everything is clean, it stays quiet.

### Step 3: Rendering & Delivery — From Markdown to Live Page

The AI writes the output as a Markdown file with YAML front matter and saves it to my [Jekyll website repo](https://github.com/CHENyiru3.github.io) under `_research_radar/YYYY-MM-DD.md`. The Hermes cron job commits and pushes this file to GitHub.

From there, GitHub Pages rebuilds the site — about 40 seconds from push to live. Jekyll's `collections` feature treats the digest directory as a content collection, and a set of Liquid templates renders each digest into the index page, individual permalink pages, and an RSS feed. Voilà.

## Design Choices

A few decisions that shaped this system:

**SQLite over API**. NetNewsWire's local database is where I track what I've read and what's new. No API rate limits, no auth tokens to rotate, no hoping a third-party service stays up. If NetNewsWire works, the pipeline works.

**AI as curator, not filter**. The model doesn't just say "relevant" or "not relevant." It generates explanations — *why* a paper matters and *why* it matters to me. That turns the digest from a sterile list into something I can actually think with. I can skip a paper and still know what I'm skipping.

**Separation of concerns**. The AI generates content. GitHub Pages renders presentation. The QA job keeps things correct. None of these layers know about each other. If I swap DeepSeek for another model, the rendering pipeline doesn't change. If I redesign the website, the AI prompt doesn't change.

**Conservative automation**. The system adds articles but never removes them. It labels, ranks, and explains — but I'm the one who decides what actually goes into my Zotero library and what shapes my research direction. Automation amplifies my attention; it doesn't replace it. And it never should.

## The Takeaway

The Research Radar isn't a complicated system — it's a cron job, an API call, and a static site generator wired together. But the *idea* is worth saying out loud: **automation that handles triage, not consumption**. The machine decides what's noise so you can decide what's signal.

If you're drowning in academic literature, the answer isn't to read faster. It's to build a better filter — tuned to your own research area, with a source list you enrich step by step. Hope this helps you get started.

  </section>

  <section id="research-radar-zh" class="bilingual-post-panel" data-bilingual-lang="zh" lang="zh-CN" hidden markdown="1">

嗨，各位博客读者！👋 欢迎来到互联网的这个角落，很高兴你路过。

每天早上 09:15，我还在喝第一杯咖啡的时候，一个我取名叫 **Clawdie** 的 AI 代理（名字来自 Claude Code，第一个带我入坑的 AI 代理）已经开始扫描二十几个期刊刚出炉的论文，读它们的摘要，然后决定哪 13 到 15 篇值得我关注。Clawdie 跑在 [Hermes Agent](https://github.com/nousresearch/hermes-agent) 上当调度大脑，cron 定时任务管排期，真正的阅读和推理丢给 **DeepSeek-V4-Pro**（爱死 DS 干这种活了）。

等我打开电脑的时候，一份排好序的简报已经在我的 [Research Radar 页面](https://chenyiru3.github.io/research-radar/)上等着我了——摘要、相关性解释，还有一个串联起当日最强主题的"热点话题"。我扫一眼列表和推荐理由，看到跟我的研究或兴趣沾边的文章，就点进 DOI 链接读原文。注意力是我唯一需要的，Research Radar 帮我省点时间。

## 为什么我要搭这个

学术出版太快了。光是 Nature、Science、Cell 及其子刊，每天就有源源不断的论文——这还没算预印本。你根本读不完。而且大部分跟我的研究方向无关。

但有时候，有那么几篇是至关重要的。关键信息被噪音淹没了，错过真的很可惜。

![我每天在 macOS 上用的 RSS 阅读器——NetNewsWire。](https://i.imgur.com/lAzG0vx.png)

我试过 RSS 阅读器（顺便说一句，我用的就是 [NetNewsWire](https://github.com/Ranchero-Software/NetNewsWire/)）、邮件订阅提醒、Zotero 关键词搜索。甚至社交媒体——Twitter、微信、小红书——但你懂的。它们全都有同一个毛病：需要*我*来筛选。每一天。天天对着 50+ 篇论文点"不相关、不相关、不相关"，只为找到那 3 篇重要的——累死个人。

Research Radar 把它翻了过来——机器做筛选，我做阅读。

## 每日简报长什么样

聊怎么搭之前，先看看每天早上页面上到底会出现什么。

**热点话题** —— Clawdie 结合我最近的阅读模式、Zotero 新增条目，加上系统级的科研方向提示，把当天 2–5 篇选中论文之间的最强主线拉出来。给我一段叙述总结，再加几条信号要点——相当于用一段话回答"今天该关注什么？"

**策展文章分区** —— 简报分成三个版块（周五多加一个）：

- **计算**（5 篇）—— 方法、算法和 AI 工具
- **生物医学**（5 篇）—— 跟我研究相关的生物学发现
- **其他领域**（3 篇）—— 生物医学之外的 AI 突破，值得瞄一眼
- **Biotech 新闻速递**（5 篇，仅周五）—— 行业动态

每篇文章我都能拿到：简洁的贡献摘要、"为什么重要"的领域解释、针对我的个性化"为什么对 Yiru 重要"，以及推荐等级（精读 / 浏览 / 知晓）。有东西引起我注意了，就点 DOI 链接去读真正的论文。

![Research Radar 示例——一份包含热点话题和策展文章分区的每日简报。](https://i.imgur.com/FfNakde.png)

## 管线：怎么复现

三步。你照着来就能搭自己的。

### 第一步：信息源收集——管好你的来源质量

来自 20+ 个期刊和预印本服务器的精选 RSS 源汇入 [NetNewsWire](https://netnews.wire.com/)，一款 macOS 原生 RSS 阅读器。我按文件夹组织这些源：Nature 系列期刊的 AOP（已接收在线论文）、Cell Press 期刊、Science，再加上 PLOS Computational Biology 和 bioRxiv 方法区这类计算生物学阵地。

NetNewsWire 把文章数据库存成本地 SQLite 文件。这为什么重要？因为整个信息源状态都可以直接查，不用 API 密钥、不用 OAuth 令牌、不依赖任何第三方服务。没有账号要管，没有 token 要轮换。就是磁盘上的一个 SQLite 数据库。

### 第二步：AI 策展——让模型来读、来排

09:15，一个 [Hermes Agent](https://hermes-agent.nousresearch.com/) cron 任务醒过来，连上 NetNewsWire 的 SQLite 数据库，抓出过去 24 小时所有未读文章，然后丢给 DeepSeek-V4-Pro，附带一份详细的策展提示词。

提示词编码了我的研究画像——目前：空间组学、单细胞分析、计算免疫学、肿瘤微环境生物学、生物医学 AI 和大模型基础模型。模型为每篇文章打相关性分，给每个版块挑最好的文章，生成你上面看到的所有摘要、解释、热点话题和推荐标签。

幕后，第二个 cron 任务在 10:00 启动，当 QA 看门狗。AI 生成的 YAML 很脆——早期有一次，缺了个闭合的 `---`，整份简报直接渲染成空白页。看门狗检查 YAML 结构（分隔符没丢、DOI 没空、数组格式对），在问题到达线上站点之前自动修掉。一切正常的时候，它不出声。

### 第三步：渲染与交付——从 Markdown 到线上页面

AI 把输出写成带 YAML 前置元数据的 Markdown 文件，存到我的 [Jekyll 网站仓库](https://github.com/CHENyiru3.github.io)里的 `_research_radar/YYYY-MM-DD.md`。Hermes cron 任务把它 commit 并 push 到 GitHub。

从那里，GitHub Pages 重建整个站点——推送到上线大概 40 秒。Jekyll 的 `collections` 把简报目录当内容集合，一组 Liquid 模板把每份简报渲染成索引页、独立永久链接页和 RSS 订阅源。就好了。

## 设计选择

几个决定了这个系统样子的选择：

**SQLite 而不是 API**。NetNewsWire 的本地数据库就是我的"读了什么、什么是新的"的唯一记录。没有 API 速率限制、不用轮换 token、不用祈祷第三方服务别挂。NetNewsWire 正常，管线就正常。

**AI 当策展人，不当过滤器**。模型不光说"相关"或"不相关"。它生成解释——*为什么*一篇论文重要，*为什么*它对我重要。这把简报从干巴巴的列表变成了我真正能拿来想事情的东西。跳过一篇论文的同时，我依然知道自己在跳过什么。

**各管各的**。AI 生成内容。GitHub Pages 渲染展示。QA 任务确保正确。这几个层次互不相知。我把 DeepSeek 换成别的模型，渲染管线不动。我重新设计网站，AI 提示词不动。

**保守地自动化**。系统只添加文章，从不删除。它贴标签、排序、解释——但什么真正进我的 Zotero 库、什么影响我的研究方向，我来定。自动化放大我的注意力，不替代它。而且永远不该替代。

## 结语

Research Radar 不是什么复杂系统——一个 cron 任务、一次 API 调用、一个静态网站生成器串在一起。但这个*想法*值得讲清楚：**自动化管的是筛选，不是阅读**。机器定什么是噪音，你来定什么是信号。

如果你被学术文献淹没了，答案不是读得更快。是搭一个更好的过滤器——调到你自己的研究领域，一步步养肥你的信息源列表。希望这篇文章对你有用！

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
  width: auto;
  max-width: min(100%, calc(100vw - 2rem));
  max-height: 72vh;
  height: auto;
  object-fit: contain;
  margin: 1.25rem auto;
  border-radius: 8px;
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
