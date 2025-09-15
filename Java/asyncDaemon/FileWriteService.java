
@Service
public class FileWriteService {
    private static final Logger logger = LoggerFactory.getLogger(FileWriteService.class);
    
    @Async("fileWriteExecutor")
    public void continuousWriteToFile(String filePath) {
        try {
            Path path = Paths.get(filePath);
            if (!Files.exists(path)) {
                Files.createFile(path);
            }
            
            AsynchronousFileChannel fileChannel = AsynchronousFileChannel.open(
                path, StandardOpenOption.WRITE);
                
            while (true) {
                String content = "Log entry: " + LocalDateTime.now() + "\n";
                ByteBuffer buffer = ByteBuffer.wrap(content.getBytes());
                
                Future<Integer> operation = fileChannel.write(buffer, fileChannel.size());
                while (!operation.isDone()) {
                    // 等待写入完成
                }
                
                buffer.clear();
                logger.info("Content written to file by {}", Thread.currentThread().getName());
                Thread.sleep(1000);
            }
        } catch (Exception e) {
            logger.error("Error in file writing task", e);
        }
    }
}
