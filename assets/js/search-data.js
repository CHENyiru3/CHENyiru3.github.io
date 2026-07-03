// get the ninja-keys element
const ninja = document.querySelector('ninja-keys');

// add the home and posts menu items
ninja.data = [{
    id: "nav-about",
    title: "About",
    section: "Navigation",
    handler: () => {
      window.location.href = "/";
    },
  },{id: "nav-projects",
          title: "Projects",
          description: "Research projects, posters, and coursework.",
          section: "Navigation",
          handler: () => {
            window.location.href = "/projects/";
          },
        },{id: "nav-research-radar",
          title: "Research Radar",
          description: "Daily academic article recommendations from curated scholarly feeds.",
          section: "Navigation",
          handler: () => {
            window.location.href = "/research-radar/";
          },
        },{id: "nav-tech-blog",
          title: "Tech Blog",
          description: "Notes on tooling, engineering, and computational workflows.",
          section: "Navigation",
          handler: () => {
            window.location.href = "/tech-blog/";
          },
        },{id: "nav-life",
          title: "Life",
          description: "Beyond the lab — interests, hobbies, and what drives me",
          section: "Navigation",
          handler: () => {
            window.location.href = "/life/";
          },
        },{id: "post-ai-engineering-mindsets-and-practices-for-computational-biologists-part-i",
        
          title: "AI Engineering Mindsets and Practices for Computational Biologists (Part I)",
        
        description: "Exploring how Spec-Driven Development shifts AI-assisted coding from vibe coding to controlled system building, and why this engineering mindset matters deeply for computational biologists.",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2026/ai-engineering-mindsets-computational-biologists-part-1/";
          
        },
      },{id: "post-how-the-research-radar-works-an-automated-literature-curation-pipeline-not-just-for-fun",
        
          title: "How the Research Radar Works: An Automated Literature Curation Pipeline (Not Just for...",
        
        description: "A behind-the-scenes look at the automated pipeline that scans academic RSS feeds, curates recommendations with AI, and renders a daily research digest — from cron jobs to Jekyll templates.",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2026/how-research-radar-works/";
          
        },
      },{id: "post-a-post-that-can-be-cited",
        
          title: "a post that can be cited",
        
        description: "this is what a post that can be cited looks like",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2024/post-citation/";
          
        },
      },{id: "post-a-post-with-tikzjax",
        
          title: "a post with TikZJax",
        
        description: "this is what included TikZ code could look like",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2023/tikzjax/";
          
        },
      },{id: "post-post-bibliography",
        
          title: "Post Bibliography",
        
        description: "",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2023/post-bibliography/";
          
        },
      },{id: "post-a-post-with-jupyter-notebook",
        
          title: "a post with jupyter notebook",
        
        description: "an example of a blog post with jupyter notebook",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2023/jupyter-notebook/";
          
        },
      },{id: "post-a-post-with-table-of-contents-on-a-sidebar",
        
          title: "a post with table of contents on a sidebar",
        
        description: "an example of a blog post with table of contents on a sidebar",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2023/sidebar-table-of-contents/";
          
        },
      },{id: "post-displaying-beautiful-tables-with-bootstrap-tables",
        
          title: "displaying beautiful tables with Bootstrap Tables",
        
        description: "an example of how to use Bootstrap Tables",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2023/tables/";
          
        },
      },{id: "post-a-post-with-table-of-contents",
        
          title: "a post with table of contents",
        
        description: "an example of a blog post with table of contents",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2023/table-of-contents/";
          
        },
      },{id: "post-a-post-with-math",
        
          title: "a post with math",
        
        description: "an example of a blog post with some math",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2015/math/";
          
        },
      },{id: "books-the-godfather",
          title: 'The Godfather',
          description: "",
          section: "Books",handler: () => {
              window.location.href = "/books/the_godfather/";
            },},{id: "news-started-a-new-independent-research-project-on-spatialtcr-focused-on-building-an-integrated-platform-for-spatial-t-cell-receptor-sequencing",
          title: '🚀 Started a new independent research project on SpatialTCR, focused on building an...',
          description: "",
          section: "News",},{id: "news-spatialmeta-paper-accepted",
          title: 'SpatialMETA Paper Accepted',
          description: "",
          section: "News",handler: () => {
              window.location.href = "/news/announcement_2/";
            },},{id: "news-best-poster-presentation-award-at-gpb-omics-amp-amp-bioinformatics-frontiers-symposium-for-our-spatialtcr-project",
          title: '🏆 Best Poster Presentation Award at GPB Omics &amp;amp;amp; Bioinformatics Frontiers Symposium for...',
          description: "",
          section: "News",},{id: "news-duke-nus-phd-offer-accepted",
          title: 'Duke-NUS PhD Offer Accepted',
          description: "",
          section: "News",handler: () => {
              window.location.href = "/news/announcement_4/";
            },},{id: "projects-spatialmeta",
          title: 'SpatialMETA',
          description: "Integrating Cross-Sample and Cross-Modality Data for Spatial Transcriptomics and Metabolomics with CVAE",
          section: "Projects",handler: () => {
              window.location.href = "/projects/1_spatialmeta/";
            },},{id: "projects-spatialtcr",
          title: 'SpatialTCR',
          description: "An integrated platform for high-resolution spatial sequencing of T cell receptor repertoires",
          section: "Projects",handler: () => {
              window.location.href = "/projects/2_spatialtcr/";
            },},{id: "projects-feast",
          title: 'FEAST',
          description: "Simulation and interpolation of spatial transcriptomics from parameter cloud",
          section: "Projects",handler: () => {
              window.location.href = "/projects/3_stmulator/";
            },},{id: "projects-temsomap",
          title: 'TemSOMap',
          description: "Mapping lineage-resolved scRNA-seq data with spatial transcriptomics",
          section: "Projects",handler: () => {
              window.location.href = "/projects/4_temsomap/";
            },},{id: "research_radar-research-radar-2026-04-28",
          title: 'Research Radar — 2026-04-28',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-04-28/";
            },},{id: "research_radar-research-radar-2026-04-29",
          title: 'Research Radar — 2026-04-29',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-04-29/";
            },},{id: "research_radar-research-radar-2026-04-30",
          title: 'Research Radar — 2026-04-30',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-04-30/";
            },},{id: "research_radar-research-radar-2026-05-01",
          title: 'Research Radar — 2026-05-01',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-01/";
            },},{id: "research_radar-research-radar-2026-05-02",
          title: 'Research Radar — 2026-05-02',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-02/";
            },},{id: "research_radar-research-radar-2026-05-03",
          title: 'Research Radar — 2026-05-03',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-03/";
            },},{id: "research_radar-research-radar-2026-05-04",
          title: 'Research Radar — 2026-05-04',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-04/";
            },},{id: "research_radar-research-radar-2026-05-05",
          title: 'Research Radar — 2026-05-05',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-05/";
            },},{id: "research_radar-research-radar-2026-05-06",
          title: 'Research Radar — 2026-05-06',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-06/";
            },},{id: "research_radar-research-radar-2026-05-07",
          title: 'Research Radar — 2026-05-07',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-07/";
            },},{id: "research_radar-research-radar-2026-05-08",
          title: 'Research Radar — 2026-05-08',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-08/";
            },},{id: "research_radar-research-radar-2026-05-09",
          title: 'Research Radar — 2026-05-09',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-09/";
            },},{id: "research_radar-research-radar-2026-05-10",
          title: 'Research Radar — 2026-05-10',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-10/";
            },},{id: "research_radar-research-radar-2026-05-11",
          title: 'Research Radar — 2026-05-11',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-11/";
            },},{id: "research_radar-research-radar-2026-05-12",
          title: 'Research Radar — 2026-05-12',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-12/";
            },},{id: "research_radar-research-radar-2026-05-14",
          title: 'Research Radar — 2026-05-14',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-14/";
            },},{id: "research_radar-research-radar-2026-05-15",
          title: 'Research Radar — 2026-05-15',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-15/";
            },},{id: "research_radar-research-radar-2026-05-16",
          title: 'Research Radar — 2026-05-16',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-16/";
            },},{id: "research_radar-research-radar-2026-05-17",
          title: 'Research Radar — 2026-05-17',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-17/";
            },},{id: "research_radar-research-radar-2026-05-18",
          title: 'Research Radar — 2026-05-18',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-18/";
            },},{id: "research_radar-research-radar-2026-05-19",
          title: 'Research Radar — 2026-05-19',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-19/";
            },},{id: "research_radar-research-radar-2026-05-20",
          title: 'Research Radar — 2026-05-20',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-20/";
            },},{id: "research_radar-research-radar-2026-05-21",
          title: 'Research Radar — 2026-05-21',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-21/";
            },},{id: "research_radar-research-radar-2026-05-23",
          title: 'Research Radar — 2026-05-23',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-23/";
            },},{id: "research_radar-research-radar-2026-05-24",
          title: 'Research Radar — 2026-05-24',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-24/";
            },},{id: "research_radar-research-radar-2026-05-25",
          title: 'Research Radar — 2026-05-25',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-25/";
            },},{id: "research_radar-research-radar-2026-05-27",
          title: 'Research Radar — 2026-05-27',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-27/";
            },},{id: "research_radar-research-radar-2026-05-28",
          title: 'Research Radar — 2026-05-28',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-28/";
            },},{id: "research_radar-research-radar-2026-05-29",
          title: 'Research Radar — 2026-05-29',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-29/";
            },},{id: "research_radar-research-radar-2026-05-31",
          title: 'Research Radar — 2026-05-31',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-05-31/";
            },},{id: "research_radar-research-radar-2026-06-01",
          title: 'Research Radar — 2026-06-01',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-01/";
            },},{id: "research_radar-research-radar-2026-06-02",
          title: 'Research Radar — 2026-06-02',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-02/";
            },},{id: "research_radar-research-radar-2026-06-03",
          title: 'Research Radar — 2026-06-03',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-03/";
            },},{id: "research_radar-research-radar-2026-06-04",
          title: 'Research Radar — 2026-06-04',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-04/";
            },},{id: "research_radar-research-radar-2026-06-07",
          title: 'Research Radar — 2026-06-07',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-07/";
            },},{id: "research_radar-research-radar-2026-06-10",
          title: 'Research Radar — 2026-06-10',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-10/";
            },},{id: "research_radar-research-radar-2026-06-11",
          title: 'Research Radar — 2026-06-11',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-11/";
            },},{id: "research_radar-research-radar-2026-06-12",
          title: 'Research Radar — 2026-06-12',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-12/";
            },},{id: "research_radar-research-radar-2026-06-13",
          title: 'Research Radar — 2026-06-13',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-13/";
            },},{id: "research_radar-research-radar-2026-06-14",
          title: 'Research Radar — 2026-06-14',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-14/";
            },},{id: "research_radar-research-radar-2026-06-15",
          title: 'Research Radar — 2026-06-15',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-15/";
            },},{id: "research_radar-research-radar-2026-06-17",
          title: 'Research Radar — 2026-06-17',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-17/";
            },},{id: "research_radar-research-radar-2026-06-18",
          title: 'Research Radar — 2026-06-18',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-18/";
            },},{id: "research_radar-research-radar-2026-06-20",
          title: 'Research Radar — 2026-06-20',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-20/";
            },},{id: "research_radar-research-radar-2026-06-21",
          title: 'Research Radar — 2026-06-21',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-21/";
            },},{id: "research_radar-research-radar-2026-06-22",
          title: 'Research Radar — 2026-06-22',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-22/";
            },},{id: "research_radar-research-radar-2026-06-24",
          title: 'Research Radar — 2026-06-24',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-24/";
            },},{id: "research_radar-research-radar-2026-06-25",
          title: 'Research Radar — 2026-06-25',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-25/";
            },},{id: "research_radar-research-radar-2026-06-26",
          title: 'Research Radar — 2026-06-26',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-26/";
            },},{id: "research_radar-research-radar-2026-06-27",
          title: 'Research Radar — 2026-06-27',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-27/";
            },},{id: "research_radar-research-radar-2026-06-29",
          title: 'Research Radar — 2026-06-29',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-29/";
            },},{id: "research_radar-research-radar-2026-06-30",
          title: 'Research Radar — 2026-06-30',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-06-30/";
            },},{id: "research_radar-research-radar-2026-07-01",
          title: 'Research Radar — 2026-07-01',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-07-01/";
            },},{id: "research_radar-research-radar-2026-07-02",
          title: 'Research Radar — 2026-07-02',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-07-02/";
            },},{id: "research_radar-research-radar-2026-07-03",
          title: 'Research Radar — 2026-07-03',
          description: "",
          section: "Research_radar",handler: () => {
              window.location.href = "/research-radar/2026-07-03/";
            },},{
        id: 'social-email',
        title: 'email',
        section: 'Socials',
        handler: () => {
          window.open("mailto:%79%69%72%75.%32%32@%69%6E%74%6C.%7A%6A%75.%65%64%75.%63%6E", "_blank");
        },
      },{
        id: 'social-github',
        title: 'GitHub',
        section: 'Socials',
        handler: () => {
          window.open("https://github.com/CHENyiru3", "_blank");
        },
      },{
        id: 'social-linkedin',
        title: 'LinkedIn',
        section: 'Socials',
        handler: () => {
          window.open("https://www.linkedin.com/in/yiru-chen-a558163b4", "_blank");
        },
      },{
        id: 'social-orcid',
        title: 'ORCID',
        section: 'Socials',
        handler: () => {
          window.open("https://orcid.org/0009-0002-5114-4947", "_blank");
        },
      },{
        id: 'social-rss',
        title: 'RSS Feed',
        section: 'Socials',
        handler: () => {
          window.open("/feed.xml", "_blank");
        },
      },{
        id: 'social-scholar',
        title: 'Google Scholar',
        section: 'Socials',
        handler: () => {
          window.open("https://scholar.google.com/citations?user=Rfv54HwAAAAJ", "_blank");
        },
      },{
      id: 'light-theme',
      title: 'Change theme to light',
      description: 'Change the theme of the site to Light',
      section: 'Theme',
      handler: () => {
        setThemeSetting("light");
      },
    },
    {
      id: 'dark-theme',
      title: 'Change theme to dark',
      description: 'Change the theme of the site to Dark',
      section: 'Theme',
      handler: () => {
        setThemeSetting("dark");
      },
    },
    {
      id: 'system-theme',
      title: 'Use system default theme',
      description: 'Change the theme of the site to System Default',
      section: 'Theme',
      handler: () => {
        setThemeSetting("system");
      },
    },];
