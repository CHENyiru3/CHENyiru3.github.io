require 'open3'

module Jekyll
  module PageLastUpdated
    module_function

    def apply(item, site)
      return if item.data['show_page_last_updated'] == false

      item.data['page_last_updated'] =
        explicit_update(item) || git_update(item, site) || file_update(item, site) || site.time
    end

    def explicit_update(item)
      item.data['page_last_updated'] || item.data['last_updated']
    end

    def git_update(item, site)
      relative_path = source_relative_path(item)
      return unless relative_path

      stdout, status = Open3.capture2(
        'git', '-C', site.source, 'log', '-1', '--format=%cs', '--', relative_path
      )
      date = stdout.strip
      status.success? && !date.empty? ? date : nil
    rescue StandardError
      nil
    end

    def file_update(item, site)
      relative_path = source_relative_path(item)
      return unless relative_path

      path = site.in_source_dir(relative_path)
      File.exist?(path) ? File.mtime(path) : nil
    end

    def source_relative_path(item)
      path = nil
      if item.respond_to?(:relative_path)
        path = item.relative_path
      end
      path = item.path if (path.nil? || path.empty?) && item.respond_to?(:path)
      path unless path.nil? || path.empty?
    end
  end
end

Jekyll::Hooks.register [:pages, :documents], :pre_render do |item, payload|
  Jekyll::PageLastUpdated.apply(item, payload['site'])
end
