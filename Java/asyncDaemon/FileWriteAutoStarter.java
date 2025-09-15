
@Component
@Order(1)
public class FileWriteAutoStarter implements ApplicationRunner {
    private static final Logger logger = LoggerFactory.getLogger(FileWriteAutoStarter.class);
    
    @Autowired
    private FileWriteService fileWriteService;
    
    @Override
    public void run(ApplicationArguments args) throws Exception {
        logger.info("Application started - auto starting file writing task via ApplicationRunner");
        fileWriteService.continuousWriteToFile("auto-start.log");
    }
}
